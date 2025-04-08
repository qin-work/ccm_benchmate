import os
import subprocess
import warnings

import torch
from Bio import Seq, SeqIO
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

from ccm_demo.sequence.utils import *


class Sequence:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.device="cuda" if torch.cuda.is_available() else "cpu"

    def embeddings(self, model="esmc_300", normalize=False):
        if model not in ["esmc_300", "esmc_600"]:
            raise ValueError("Invalid model name")
        protein=ESMProtein(self.sequence)
        client = ESMC.from_pretrained(model).to(self.device)  # or "cpu"
        protein_tensor = client.encode(protein)
        logits_output = client.logits(
            protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
        )
        if normalize:
            logits_output=logits_output[0].mean(dim=0)
        return logits_output

    def mutate(self, position, to, inplace=False, new_name=name):
        """
        position is 0 based
        """
        new_sequence = self.sequence
        new_sequence[position] = to

        if not inplace:
            new_instance = Sequence(new_name, new_sequence, self.sequence)
            return new_instance
        else:
            self.sequence = new_sequence
            self.name = new_name
            return self

    def msa(self, database, destination):
        if not os.path.exists(database):
            raise FileNotFoundError("There is no such database: {}".format(database))

        self.write(os.path.join(os.getcwd(), "query.fasta"))

        run=subprocess.run(["run_mmseqs.sh" " -i", str(os.path.join(os.getcwd(), "query.fasta")),
                            "-d", database, "-o", destination, "-c"])
        if run.returncode != 0:
            warnings.warn("There was an error with mmseqs run see error below '\n' {}".format(run.stderr))
        else:
            print("mmseqs run complese the results are in {}".format(destination))

    #TODO Biotite
    def align(self, other):
        pass

    def write(self, path):
        seq=SeqIO.SeqRecord(Seq(self.sequence), id=self.name, description="")
        with open(path, "w") as handle:
            SeqIO.write(seq, handle, "fasta")

    def __str__(self):
        print("Sequence with name {} and {} aas".format(self.name, len(self.sequence)))

