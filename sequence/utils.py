import logging
import os
import random
import tarfile
import time
from typing import Union

import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)
host_url="https://api.colabfold.com"

headers = {}
headers["User-Agent"] = "ccm"

TQDM_BAR_FORMAT = (
    "{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]"
)

def submit(seqs, use_pairing=None, use_filter=True, pairing_strategy="greedy", use_env=True, N=101,
           host_url=host_url, headers=headers):
    """
    :param seqs: sequences
    :param mode: paired or single
    :param N: max number sequences to return
    :return: submission status and ticket id if successful
    """
    n, query = N, ""
    for seq in seqs:
        query += f">{n}\n{seq}\n"
        n += 1


    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    if len(seqs) < 1 and use_pairing is None:
        use_pairing = True
    if len(seqs) ==1:
        use_pairing = False

    if use_pairing:
        mode = ""
        # greedy is default, complete was the previous behavior
        if pairing_strategy == "greedy":
            mode = "pairgreedy"
        elif pairing_strategy == "complete":
            mode = "paircomplete"
        if use_env:
            mode = mode + "-env"

    submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"

    while True:
        error_count = 0
        try:
            res = requests.post(
                f"{host_url}/{submission_endpoint}",
                data={"q": query, "mode": mode},
                timeout=6.02,
                headers=headers,
            )
        except requests.exceptions.Timeout:
            logger.warning("Timeout while submitting to MSA server. Retrying...")
            continue
        except Exception as e:
            error_count += 1
            logger.warning(
                f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
            )
            logger.warning(f"Error: {e}")
            time.sleep(5)
            if error_count > 5:
                raise
            continue
        break

    try:
        out = res.json()
    except ValueError:
        logger.error(f"Server didn't reply with json: {res.text}")
        out = {"status": "ERROR"}
    return out

def status(ID, host_url=host_url, headers=headers):
    while True:
        error_count = 0
        try:
            res = requests.get(
                f"{host_url}/ticket/{ID}", timeout=6.02, headers=headers
            )
        except requests.exceptions.Timeout:
            logger.warning(
                "Timeout while fetching status from MSA server. Retrying..."
            )
            continue
        except Exception as e:
            error_count += 1
            logger.warning(
                f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
            )
            logger.warning(f"Error: {e}")
            time.sleep(5)
            if error_count > 5:
                raise
            continue
        break
    try:
        out = res.json()
    except ValueError:
        logger.error(f"Server didn't reply with json: {res.text}")
        out = {"status": "ERROR"}
    return out

def download(ID, path, host_url=host_url, headers=headers):
    error_count = 0

    if not path.endswith(".tar.gz"):
        path = f"{path}.tar.gz"

    while True:
        try:
            res = requests.get(
                f"{host_url}/result/download/{ID}", timeout=6.02, headers=headers
            )
        except requests.exceptions.Timeout:
            logger.warning(
                "Timeout while fetching result from MSA server. Retrying..."
            )
            continue
        except Exception as e:
            error_count += 1
            logger.warning(
                f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
            )
            logger.warning(f"Error: {e}")
            time.sleep(5)
            if error_count > 5:
                raise
            continue
        break
    with open(path, "wb") as out:
        out.write(res.content)

def process(sequences, path, use_env=True, use_filter=True, use_pairing=None, pairing_strategy="greedy"):
    # process input x
    seqs = [sequences] if isinstance(sequences, str) else sequences

    # setup mode
    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    if use_pairing:
        mode = ""
        # greedy is default, complete was the previous behavior
        if pairing_strategy == "greedy":
            mode = "pairgreedy"
        elif pairing_strategy == "complete":
            mode = "paircomplete"
        if use_env:
            mode = mode + "-env"

    # call mmseqs2 api
    tar_gz_file = path
    N, REDO = 101, True

    # deduplicate and keep track of order
    seqs_unique = []
    # TODO this might be slow for large sets
    [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
    Ms = [N + seqs_unique.index(seq) for seq in seqs]
    # lets do it!
    if not os.path.isfile(tar_gz_file):
        TIME_ESTIMATE = 150 * len(seqs_unique)
        with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
            while REDO:
                pbar.set_description("SUBMIT")

                # Resubmit job until it goes through
                out = submit(seqs_unique, mode, N)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                    sleep_time = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
                    # resubmit
                    time.sleep(sleep_time)
                    out = submit(seqs_unique, mode, N)

                if out["status"] == "ERROR":
                    msg = (
                        "MMseqs2 API is giving errors. Please confirm your "
                        " input is a valid protein sequence. If error persists, "
                        "please try again an hour later."
                    )
                    raise Exception(msg)

                if out["status"] == "MAINTENANCE":
                    msg = (
                        "MMseqs2 API is undergoing maintenance. "
                        "Please try again in a few minutes."
                    )
                    raise Exception(msg)

                # wait for job to finish
                ID, TIME = out["id"], 0
                pbar.set_description(out["status"])
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                    t = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
                    time.sleep(t)
                    out = status(ID)
                    pbar.set_description(out["status"])
                    if out["status"] == "RUNNING":
                        TIME += t
                        pbar.update(n=t)

                if out["status"] == "COMPLETE":
                    if TIME < TIME_ESTIMATE:
                        pbar.update(n=(TIME_ESTIMATE - TIME))
                    REDO = False

                if out["status"] == "ERROR":
                    REDO = False
                    msg = (
                        "MMseqs2 API is giving errors. Please confirm your "
                        " input is a valid protein sequence. If error persists, "
                        "please try again an hour later."
                    )
                    raise Exception(msg)

            # Download results
            download(ID, tar_gz_file)

    # prep list of a3m files
    if use_pairing:
        a3m_files = [f"{path}/pair.a3m"]
    else:
        a3m_files = [f"{path}/uniref.a3m"]
        if use_env:
            a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

    # extract a3m files
    if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(path)

    # gather a3m lines
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        for line in open(a3m_file, "r"):
            if len(line) > 0:
                if "\x00" in line:
                    line = line.replace("\x00", "")
                    update_M = True
                if line.startswith(">") and update_M:
                    M = int(line[1:].rstrip())
                    update_M = False
                    if M not in a3m_lines:
                        a3m_lines[M] = []
                a3m_lines[M].append(line)

    a3m_lines = ["".join(a3m_lines[n]) for n in Ms]
    return a3m_lines



