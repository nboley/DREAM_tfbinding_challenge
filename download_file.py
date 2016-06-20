""" Example file submission script

Requires the `aws` command line utility: http://aws.amazon.com/cli/
"""
import hashlib
import json
import os
import requests
import subprocess
import sys
import time
import urlparse

from pyTFbindtools.ENCODE_ChIPseq_tools import (
    find_ENCODE_DCC_experiment_files_from_experiment_json,
    extract_ENCODE_DCC_fastqs)

experiment_urls = """
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR014YCR
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR047BUZ
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR545FXC
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR267DFA
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR066EBK
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR061DGF
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR869IUD
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR857SDB
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR563LLO
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR189TRZ
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR888XZK
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR736BUG
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR290ZOS
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR735KEY
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR324RCI
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR310NYI
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR080XEY
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR038GMB
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR350ORK
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR601OGE
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR445QRF
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR521IID
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR893QWP
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR867WPH
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR062QVQ
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR016BMM
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR881UOO
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR866QPZ
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR897JAS
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR569XNP
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR126YEB
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR359LOD
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR932ZMX
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR107PGG
""".split()

# these are the adult liver controls
experiment_urls = """
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR621CSR/
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR019XRC/
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR628ONI/
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR931XLP/
""".split()

# these are the child liver controls
experiment_urls = """
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR443AIH/
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR769JUI/
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR849OTN/
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR904IXG/
""".split()

# IPSC and 
experiment_urls = """
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR447FVP/
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR655HDE/
https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/ENCSR328AHC/
""".split()

experiment_urls = """
https://www.encodeproject.org/experiments/ENCSR195ZCD/
https://www.encodeproject.org/experiments/ENCSR398JTO/
""".split()

encoded_access_key = 'ONRL3ZRO'
encoded_secret_access_key = 'sw5vjtgc7k2fen7x'
#encoded_access_key = '2KMYOYKL'
#encoded_secret_access_key = 'pynx26slwsehxyln'

def iter_fastq_fnames(chipseq_fastqs):
    assert len(chipseq_fastqs[1]) == len(chipseq_fastqs[2])
    for x in chipseq_fastqs[0]:
        ofname = (
            "CHIPseq.%s.%s.EXPID_{0.exp_id}.BSID_{0.bsid}.BSREP{0.rep_key[0]}.TECHREP{0.rep_key[1]}.FILEID{0.id}.unpaired.fastq.gz" % (
                x.target_id.split("/")[-2].replace("-human", ""),
                x.sample_type.replace(" ", "_"))
        ).format(x)
        yield ofname, x.file_loc
    for x in chipseq_fastqs[1]:
        ofname = (
            "CHIPseq.%s.%s.EXPID_{0.exp_id}.BSID_{0.bsid}.BSREP{0.rep_key[0]}.TECHREP{0.rep_key[1]}.FILEID{0.id}.R1.fastq.gz" % (
                x.target_id.split("/")[-2].replace("-human", ""),
                x.sample_type.replace(" ", "_"))
        ).format(x)
        yield ofname, x.file_loc
    for x in chipseq_fastqs[2]:
        ofname = (
            "CHIPseq.%s.%s.EXPID_{0.exp_id}.BSID_{0.bsid}.BSREP{0.rep_key[0]}.TECHREP{0.rep_key[1]}.FILEID{0.id}.R2.fastq.gz" % (
                x.target_id.split("/")[-2].replace("-human", ""),
                x.sample_type.replace(" ", "_"))
        ).format(x)
        yield ofname, x.file_loc
    return

def download_file( url ):
    base = urlparse.urlunsplit(list(urlparse.urlsplit(url)[:2]) + ['', '', ''])
    headers = {
        'Content-type': 'application/json',
        'Accept': 'application/json',
    }

    def download_fastq(local_filename, url):
        r = requests.get(
            url,
            stream=True,
            auth=(encoded_access_key, encoded_secret_access_key),
        )
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=4096):
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
        return

    r = requests.get(
        url,
        auth=(encoded_access_key, encoded_secret_access_key),
        headers=headers,
    )
    try:
        r.raise_for_status()
    except:
        print('Request failed: %s %s' % (r.status_code, r.reason))
        if r.status_code == 409:
            print(r.text)
            return
        raise
    for ofname, url_suffix in iter_fastq_fnames(
            extract_ENCODE_DCC_fastqs(
                find_ENCODE_DCC_experiment_files_from_experiment_json(r.json())
            )):
        print "Downloading {} from {}".format(ofname, base + url_suffix)
        download_fastq(ofname, base + url_suffix)
    return

def find_controls(exp_ids):
    control_exp_ids = set()

    for experiment_id in exp_ids:
        print experiment_id
        print control_exp_ids
        #URL = "https://4024-idr203-ef5458d-jseth.demo.encodedcc.org/experiments/{}/".format(experiment_id)
        #URL = "https://www.encodeproject.org/experiments/{}/".format(experiment_id)
        response = requests.get(
            URL,
            headers={'accept': 'application/json'},
            auth=(encoded_access_key, encoded_secret_access_key)
        )
        response_json_dict = response.json()
        print response_json_dict
        control_exp_ids.update(
            x['accession'] for x in response_json_dict['possible_controls'])

def download_fastqs():
    for i, url in enumerate(experiment_urls):
        #if i != 28: continue
        print "Processing ", i, len(experiment_urls), url
        download_file(url)
        
def main():
    #find_controls([x.split("/")[-1] for x in experiment_urls])
    download_fastqs()
    
if __name__ == '__main__':
    main()
