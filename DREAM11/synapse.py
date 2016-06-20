import datetime

import urllib
import json

import synapseclient
from synapseclient import Project, Folder, File, Link, Evaluation, Wiki, Team
import synapseclient.utils as utils

syn = synapseclient.Synapse()
syn.login('nboley', 'Sp33d427')
project = syn.get('syn6044646')
wiki = syn.getWiki("syn6044646")

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
    dict(name='spearman_correlation',
         display_name='Spearman Correlation',
         columnType='DOUBLE'),
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
    return
    
def score_regression_challenge():
    pass


def main():
    #create_evaluation_queues_and_leaderboard()
    test()

if __name__ == '__main__':
    main()
