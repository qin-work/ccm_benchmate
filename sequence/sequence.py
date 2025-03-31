import copy

import torch

from Bio import Seq, SeqIO
from Bio import Blast as blast

from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig
from ccm_demo.sequence import utils


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

    def mutate(self, position, to, inplace=False, new_name=self.name):
        """
        position is 0 based
        """
        new_sequence = copy.deepcopy(self.sequence)
        new_sequence[position] = to

        if not inplace:
            new_instance = Sequence(new_name, new_sequence, self.sequence)
            return new_instance
        else:
            self.sequence = new_sequence
            self.name = new_name
            return self

    def msa(self, others=None):
        if others is None:


    def write(self):
        pass

    def __str__(self):
        print("Sequence with name {} and {} aas".format(self.name, len(self.sequence)))

