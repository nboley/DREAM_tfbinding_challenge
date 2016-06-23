import sys
import traceback

import datetime

import urllib
import json

from StringIO import StringIO

import synapseclient
from synapseclient import Project, Folder, File, Link, Evaluation, Wiki, Team
import synapseclient.utils as utils

syn = synapseclient.Synapse()
user = syn.login('nboley', 'Sp33d427')
ADMIN_USER_IDS = ["3341836",]

project = syn.get('syn6044646')
wiki = syn.getWiki("syn6044646")

## Messages for challenge scoring script.

import string
import sys
import warnings


## Module level state. You'll need to set a synapse object at least
## before using this module.
send_messages = True
send_notifications = True
acknowledge_receipt = False
dry_run = False


## Edit these URLs to point to your challenge and its support forum
defaults = dict(
    challenge_instructions_url = "https://www.synapse.org/",
    support_forum_url = "http://support.sagebase.org/sagebase")

##---------------------------------------------------------
## Message templates:
## Edit to fit your challenge.
##---------------------------------------------------------

validation_failed_subject_template = "Validation error in submission to {queue_name}"
validation_failed_template = """\
Hello {username},

Sorry, but we were unable to validate your submission to the {queue_name}.

Please refer to the challenge instructions which can be found at \
{challenge_instructions_url} and to the error message below:

submission name: {submission_name}
submission ID: {submission_id}

{message}

If you have questions, please ask on the forums at {support_forum_url}.

Sincerely,

the scoring script
"""

validation_passed_subject_template = "Submission received to {queue_name}"
validation_passed_template = """\
Hello {username},

We have received your submission to the {queue_name} and confirmed that it is correctly formatted.

submission name: {submission_name}
submission ID: {submission_id}

If you have questions, please ask on the forums at {support_forum_url} or refer to the challenge \
instructions which can be found at {challenge_instructions_url}.

Sincerely,

the scoring script
"""

scoring_succeeded_subject_template = "Scored submission to {queue_name}"
scoring_succeeded_template = """\
Hello {username},

Your submission \"{submission_name}\" (ID: {submission_id}) to the {queue_name} has been scored:

{message}

If you have questions, please ask on the forums at {support_forum_url}.

Sincerely,

the scoring script
"""

scoring_error_subject_template = "Exception while scoring submission to {queue_name}"
scoring_error_template = """\
Hello {username},

Sorry, but we were unable to process your submission to the {queue_name}.

Please refer to the challenge instructions which can be found at \
{challenge_instructions_url} and to the error message below:

submission name: {submission_name}
submission ID: {submission_id}

{message}

If you have questions, please ask on the forums at {support_forum_url}.

Sincerely,

the scoring script
"""

notification_subject_template = "Exception while scoring submission to {queue_name}"
error_notification_template = """\
Hello Challenge Administrator,

The scoring script for {queue_name} encountered an error:

{message}

Sincerely,

the scoring script
"""


class DefaultingFormatter(string.Formatter):
    """
    Python's string.format has the annoying habit of raising a KeyError
    if you don't completely fill in the template. Let's do something a
    bit nicer.
    Adapted from: http://stackoverflow.com/a/19800610/199166
    """
    def get_value(self, key, args, kwds):
        if isinstance(key, str):
            value = kwds.get(key, defaults.get(key, None))
            if value is None:
                value = "{{{0}}}".format(key)
                warnings.warn("Missing template variable %s" % value)
            return value
        else:
            Formatter.get_value(key, args, kwds)

formatter = DefaultingFormatter()

##---------------------------------------------------------
## functions for sending various types of messages
##---------------------------------------------------------

def validation_failed(userIds, **kwargs):
    if send_messages:
        return send_message(userIds=userIds, 
                            subject_template=validation_failed_subject_template,
                            message_template=validation_failed_template,
                            kwargs=kwargs)

def validation_passed(userIds, **kwargs):
    if acknowledge_receipt:
        return send_message(userIds=userIds,
                            subject_template=validation_passed_subject_template,
                            message_template=validation_passed_template,
                            kwargs=kwargs)

def scoring_succeeded(userIds, **kwargs):
    if send_messages:
        return send_message(userIds=userIds,
                            subject_template=scoring_succeeded_subject_template,
                            message_template=scoring_succeeded_template,
                            kwargs=kwargs)

def scoring_error(userIds, **kwargs):
    if send_messages:
        return send_message(userIds=userIds,
                            subject_template=scoring_error_subject_template,
                            message_template=scoring_error_template,
                            kwargs=kwargs)

def error_notification(userIds, **kwargs):
    if send_notifications:
        return send_message(userIds=userIds,
                            subject_template=notification_subject_template,
                            message_template=error_notification_template,
                            kwargs=kwargs)

def send_message(userIds, subject_template, message_template, kwargs):
    print kwargs
    subject = formatter.format(subject_template, **kwargs)
    message = formatter.format(message_template, **kwargs)
    if dry_run:
        print "\nDry Run: would have sent:"
        print subject
        print "-" * 60
        print message
        return None
    elif syn:
        response = syn.sendMessage(
            userIds=userIds,
            messageSubject=subject,
            messageBody=message)
        print "sent: ", unicode(response).encode('utf-8')
        return response
    else:
        sys.stderr.write("Can't send message. No Synapse object configured\n")
        
def get_user_name(profile):
    names = []
    if 'firstName' in profile and profile['firstName'] and profile['firstName'].strip():
        names.append(profile['firstName'])
    if 'lastName' in profile and profile['lastName'] and profile['lastName'].strip():
        names.append(profile['lastName'])
    if len(names)==0:
        names.append(profile['userName'])
    return " ".join(names)

BASE_LEADERBOARD_COLUMNS = [
    dict(name='objectId',      display_name='ID',      columnType='STRING', maximumSize=20),
    dict(name='userId',        display_name='User',    columnType='STRING', maximumSize=20, renderer='userid'),
    dict(name='entityId',      display_name='Entity',  columnType='STRING', maximumSize=20, renderer='synapseid'),
    dict(name='versionNumber', display_name='Version', columnType='INTEGER'),
    dict(name='name',          display_name='Name',    columnType='STRING', maximumSize=240),
    dict(name='team',          display_name='Team',    columnType='STRING', maximumSize=240)
]

classification_leaderboard_columns = BASE_LEADERBOARD_COLUMNS + [
    dict(name='auroc',         display_name='auROC',             columnType='DOUBLE'),
    dict(name='auprg',         display_name='auPRG',             columnType='DOUBLE'),
    dict(name='re_0_05_fdr',   display_name='Recall @ 5% FDR',   columnType='DOUBLE'),
    dict(name='re_0_10_fdr',   display_name='Recall @10% FDR',   columnType='DOUBLE'),
    dict(name='re_0_25_fdr',   display_name='Recall @25% FDR',   columnType='DOUBLE'),
    dict(name='re_0_50_fdr',   display_name='Recall @50% FDR',   columnType='DOUBLE')
]

regression_leaderboard_columns = BASE_LEADERBOARD_COLUMNS + [
    dict(name='spearman_correlation', display_name='Spearman Correlation', columnType='DOUBLE'),
    dict(name='mse', display_name='MSE', columnType='DOUBLE')
]

def create_supertable_leaderboard(evaluation, leaderboard_columns):
    """
    Create the leaderboard using a supertable, a markdown extension that dynamically
    builds a table by querying submissions. Because the supertable re-queries whenever
    the page is rendered, this step only has to be done once.
    """
    uri_base = urllib.quote_plus("/evaluation/submission/query")
    # it's incredibly picky that the equals sign here has to be urlencoded, but
    # the later equals signs CAN'T be urlencoded.
    query = urllib.quote_plus(
        'query=select * from evaluation_%s where status=="SCORED"'
        % utils.id_of(evaluation))
    params = [  ('paging', 'true'),
                ('queryTableResults', 'true'),
                ('showIfLoggedInOnly', 'false'),
                ('pageSize', '25'),
                ('showRowNumber', 'false'),
                ('jsonResultsKeyName', 'rows')]

    # Columns specifications have 4 fields: renderer, display name, column name, sort
    # Renderer and sort are usually 'none' and 'NONE'.
    for i, column in enumerate(leaderboard_columns):
        fields = {'renderer':'none', 'sort':'NONE'}
        fields.update(column)
        if 'display_name' not in fields:
            fields['display_name'] = fields['name']
        params.append((
            'columnConfig%s' % i,
            "{renderer},{display_name},{name};,{sort}".format(**fields)
        ))

    return "${supertable?path=" + uri_base + "%3F" + query + "&" + "&".join(
        [key+"="+urllib.quote_plus(value) for key,value in params]) + "}"

def delete_evaluations():
    for evaluation_name in ["Within Cell Type Classification",
                            "Between Cell Type Classification",
                            "Between Cell Type Regression"]:
        try:
            evaluation = syn.getEvaluationByName(evaluation_name)
            syn.delete(evaluation)
        except:
            pass

def create_leaderboard_wiki(
        intra_cell_type_evaluation,
        inter_cell_type_evaluation,
        regression_evaluation):
    LEADERBOARD_MARKDOWN = """\
## Within Cell Type Classification

{intra_cell_type_classification_supertable}

## Between Cell Type Classification

{inter_cell_type_classification_supertable}

## Between Cell Type Regression

{inter_cell_type_regression_supertable}
"""
    leaderboard_wiki_text = LEADERBOARD_MARKDOWN.format(
        intra_cell_type_classification_supertable=create_supertable_leaderboard(
            intra_cell_type_evaluation, classification_leaderboard_columns),
        inter_cell_type_classification_supertable=create_supertable_leaderboard(
            inter_cell_type_evaluation, classification_leaderboard_columns),
        inter_cell_type_regression_supertable=create_supertable_leaderboard(
            regression_evaluation, regression_leaderboard_columns)
    )
    lb_wiki = Wiki(
        title="Leaderboard",
        owner=project,
        parentWikiId=wiki.id,
        markdown=leaderboard_wiki_text
    )
    lb_wiki = syn.store(lb_wiki)

def find_or_create_evaluations():    
    # create the within cell type classification evaluation object
    try:
        intra_cell_type_evaluation = syn.getEvaluationByName(
            "Within Cell Type Classification")
    except:
        intra_cell_type_evaluation = syn.store(
            Evaluation(
                name="Within Cell Type Classification",
                contentSource=project.id,
                status="OPEN",
                submissionInstructionsMessage=\
                "Submit a tsv file containing predicted bin labels.",
                submissionReceiptMessage=\
                "Your submission has been received."),

            quota=dict(numberOfRounds=1,
                       roundDurationMillis=1000*60*60*48, ## 48 hours
                       submissionLimit=20,
                       firstRoundStart=datetime.datetime.now().strftime(
                           synapseclient.utils.ISO_FORMAT))
        )

    # create the within cell type classification evaluation object
    try:
        inter_cell_type_evaluation = syn.getEvaluationByName(
            "Between Cell Type Classification")
    except:
        inter_cell_type_evaluation = syn.store(
            Evaluation(
                name="Between Cell Type Classification",
                contentSource=project.id,
                status="OPEN",
                submissionInstructionsMessage=\
                "Submit a tsv file containing predicted bin labels.",
                submissionReceiptMessage=\
                "Your submission has been received."),

            quota=dict(numberOfRounds=1,
                       roundDurationMillis=1000*60*60*48, ## 48 hours
                       submissionLimit=20,
                       firstRoundStart=datetime.datetime.now().strftime(synapseclient.utils.ISO_FORMAT))
        )

    # create the within cell type classification evaluation object
    try:
        regression_evaluation = syn.getEvaluationByName(
            "Between Cell Type Regression")
    except:
        regression_evaluation = syn.store(
            Evaluation(
                name="Between Cell Type Regression",
                contentSource=project.id,
                status="OPEN",
                submissionInstructionsMessage=\
                "Submit a tsv file containing predicted average ChIP-seq read coverage.",
                submissionReceiptMessage=\
                "Your submission has been received."),

            quota=dict(numberOfRounds=1,
                       roundDurationMillis=1000*60*60*48, ## 48 hours
                       submissionLimit=20,
                       firstRoundStart=datetime.datetime.now().strftime(
                           synapseclient.utils.ISO_FORMAT))
        )
        
    return (intra_cell_type_evaluation,
            inter_cell_type_evaluation,
            regression_evaluation)

def create_evaluation_queues_and_leaderboard():
    # create evaluations
    delete_evaluations()
    evaluations = find_or_create_evaluations()
    create_leaderboard_wiki(*evaluations)

def validate_submission(evaluation, submission):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    return True, "Looks OK to me!"


def score_classification_submission(evaluation, submission):
    """
    Find the right scoring function and score the submission

    :returns: (score, message) where score is a dict of stats and message
              is text for display to user
    """
    import random
    return (dict(auroc=random.random(),
                 auprg=random.random(),
                 re_0_05_fdr=random.random(),
                 re_0_10_fdr=random.random(),
                 re_0_25_fdr=random.random(),
                 re_0_50_fdr=random.random()),
            "You did fine!")

def score_regression_submission(evaluation, submission):
    """
    Find the right scoring function and score the submission

    :returns: (score, message) where score is a dict of stats and message
              is text for display to user
    """
    import random
    return (dict(spearman_correlation=random.random(),
                 mse=random.random()),
            "You did fine!")

def validate(evaluation, dry_run=False):

    if type(evaluation) != Evaluation:
        evaluation = syn.getEvaluation(evaluation)

    print "\n\nValidating", evaluation.id, evaluation.name
    print "-" * 60
    sys.stdout.flush()

    for submission, status in syn.getSubmissionBundles(evaluation, status='RECEIVED'):
        ## refetch the submission so that we get the file path
        ## to be later replaced by a "downloadFiles" flag on getSubmissionBundles
        submission = syn.getSubmission(submission)

        print "validating", submission.id, submission.name
        try:
            is_valid, validation_message = validate_submission(evaluation, submission)
        except Exception as ex1:
            is_valid = False
            print "Exception during validation:", type(ex1), ex1, ex1.message
            traceback.print_exc()
            validation_message = str(ex1)

        status.status = "VALIDATED" if is_valid else "INVALID"

        if not dry_run:
            status = syn.store(status)

        ## send message AFTER storing status to ensure we don't get repeat messages
        profile = syn.getUserProfile(submission.userId)
        if is_valid:
            validation_passed(
                userIds=[submission.userId],
                username=get_user_name(profile),
                queue_name=evaluation.name,
                submission_id=submission.id,
                submission_name=submission.name)
        else:
            validation_failed(
                userIds=[submission.userId],
                username=get_user_name(profile),
                queue_name=evaluation.name,
                submission_id=submission.id,
                submission_name=submission.name,
                message=validation_message)


def score(evaluation, score_submission, dry_run=False):
    if type(evaluation) != Evaluation:
        evaluation = syn.getEvaluation(evaluation)
    print '\n\nScoring ', evaluation.id, evaluation.name
    print "-" * 60
    sys.stdout.flush()

    for submission, status in syn.getSubmissionBundles(evaluation, status='VALIDATED'):

        status.status = "INVALID"

        ## refetch the submission so that we get the file path
        ## to be later replaced by a "downloadFiles" flag on getSubmissionBundles
        submission = syn.getSubmission(submission)

        try:
            score, message = score_submission(evaluation, submission)

            print "scored:", submission.id, submission.name, submission.userId, score

            ## fill in team in submission status annotations
            if 'teamId' in submission:
                team = syn.restGET('/team/{id}'.format(id=submission.teamId))
                if 'name' in team:
                    score['team'] = team['name']
                else:
                    score['team'] = submission.teamId
            elif 'userId' in submission:
                profile = syn.getUserProfile(submission.userId)
                score['team'] = get_user_name(profile)
            else:
                score['team'] = '?'

            status.annotations = synapseclient.annotations.to_submission_status_annotations(score,is_private=False)
            status.status = "SCORED"
            ## if there's a table configured, update it
            #if not dry_run and evaluation.id in conf.leaderboard_tables:
            #    update_leaderboard_table(conf.leaderboard_tables[evaluation.id], submission, fields=score, dry_run=False)

        except Exception as ex1:
            sys.stderr.write('\n\nError scoring submission %s %s:\n' % (submission.name, submission.id))
            st = StringIO()
            traceback.print_exc(file=st)
            sys.stderr.write(st.getvalue())
            sys.stderr.write('\n')
            message = st.getvalue()

            if ADMIN_USER_IDS:
                submission_info = "submission id: %s\nsubmission name: %s\nsubmitted by user id: %s\n\n" % (submission.id, submission.name, submission.userId)
                error_notification(userIds=ADMIN_USER_IDS, message=submission_info+st.getvalue())

        if not dry_run:
            status = syn.store(status)

        ## send message AFTER storing status to ensure we don't get repeat messages
        profile = syn.getUserProfile(submission.userId)

        if status.status == 'SCORED':
            scoring_succeeded(
                userIds=[submission.userId],
                message=message,
                username=get_user_name(profile),
                queue_name=evaluation.name,
                submission_name=submission.name,
                submission_id=submission.id)
        else:
            scoring_failed(
                userIds=[submission.userId],
                message=message,
                username=get_user_name(profile),
                queue_name=evaluation.name,
                submission_name=submission.name,
                submission_id=submission.id)

    sys.stdout.write('\n')

    
def test():
    ids = {}
    
    # teams are really just permissions groups, and we need to create a group
    # that has permissions to submit to the challenge
    try:
        participants_team = syn.getTeam("Participants Group")
        print "Loaded team %s %s"%(participants_team.id, participants_team.name)
    except:
        participants_team = syn.store(Team(
            name='Participants Group',
            description='A team for people who have joined the challenge'))
        print syn.setPermissions(project, participants_team.id, ['READ', 'DOWNLOAD'])
        for evaluation in find_or_create_evaluations():
            print syn.setPermissions(
                evaluation,
                participants_team.id, [
                    'CREATE', 'READ', 'UPDATE',
                    'PARTICIPATE', 'SUBMIT', 'READ_PRIVATE_SUBMISSION'
                ]
            )
        print "Created team %s %s" % (participants_team.id, participants_team.name)

    ids['participants_team'] = participants_team.id

    ## the challenge object associates the challenge project with the
    ## participants team
    try:
        challenge_object = syn.restGET(
            "/entity/{ID}/challenge".format(ID=urllib.quote(project.id))
        )
        print "Loaded challenge object: %s" % challenge_object['id']
    except:
        raise
        challenge_json = {
            'participantTeamId':participants_team.id, 'projectId':project.id}
        challenge_object = syn.restPOST(
            "/challenge", body=json.dumps(challenge_json))
        print "Created challenge object: %s" % challenge_object['id']
    ids['challenge_object'] = challenge_object['id']
    
    # get a submission team
    test_submissions_teams = syn.restGET(
        '/challenge/{challengeId}/challengeTeam'.format(
            challengeId=challenge_object['id'])
    )
    if test_submissions_teams['totalNumberOfResults'] == 0:
        # create a team that will make submissions
        test_submissions_team = syn.store(Team(
            name="Test Submission Team",
            description='A team to test the submissions machinery'))
        # register team with challenge
        request_body = {
            'teamId': test_submissions_team['teamId'],
            'challengeId': challenge_object['id'] }
        test_submissions_team = syn.restPOST('/challenge/{challengeId}/challengeTeam'.format(
            challengeId=challenge_object['id']), json.dumps(request_body))
    elif test_submissions_teams['totalNumberOfResults'] >= 1:
        test_submissions_team = test_submissions_teams['results'][0]
        ids['test_submissions_team'] = test_submissions_team['id']
        print "Restored a team to test submissions: %s" % (
            test_submissions_team['id'], )
    else:
        print test_submissions_teams
        assert False, "We only expect a single test submission user."

    # Create the participant project
    try:
        participant_project = syn.get("syn6170604")
    except:
        raise
        participant_project = syn.store(
            Project(name="Test Submission Team Project"))
        print "Created project %s %s" % (
            participant_project.id, participant_project.name)

    for evaluation in find_or_create_evaluations():
        participant_file = syn.store(File(
            synapseclient.utils.make_bogus_data_file(),
            parent=participant_project))
        syn.submit(evaluation=evaluation,
                   entity=participant_file,
                   team="Participants Group",
        )
    
    for evaluation in find_or_create_evaluations():
        if evaluation.name == 'Between Cell Type Regression':
            score_submission = score_regression_submission
        else:
            score_submission = score_classification_submission
        validate(evaluation)
        score(evaluation, score_submission)
        
    return
    
def main():
    #create_evaluation_queues_and_leaderboard()
    test()

if __name__ == '__main__':
    main()
