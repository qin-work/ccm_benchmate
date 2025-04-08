import os
import subprocess

from biotite.structure import io
import torch
import requests

from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, SamplingConfig


class Structure:
    def __init__(self, pdb):
        self.pdb = pdb
        self.sasa=None,
        self.embeddings=None
        self.foldseek=None
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to(self.device)

    def download(self, id, source="PDB", destination=None, load_after_download=True):
        if source == "PDB":
            url="http://files.rcsb.org/download/{}.pdb".format(id)
        elif source == "AFDB":
            url="https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb".format(id)
        else:
            raise NotImplementedError("We can only download structures from PDB or AFDB")

        download = requests.get(url, stream=True)
        download.raise_for_status()
        with open("{}/{}.pdb".format(destination, id), "wb") as f:
            f.write(download.content)

        self.pdb="{}/{}.pdb".format(destination, id)

        if load_after_download:
            atom_array = io.load_structure(self.pdb)
        self.structure=atom_array
        return self

    def embeddings(self, model="esm3", normalize=False):
        protein = ESMProtein.from_pdb("./1utn.pdb")
        protein_tensor = self.ESM3InferenceClient.encode(protein)

        if normalize:
            output = self.ESM3InferenceClient.forward_and_sample(
                protein_tensor, SamplingConfig(return_mean_embedding=True)
            )
        else:
            output = self.ESM3InferenceClient.forward_and_sample(
                protein_tensor, SamplingConfig(return_per_residue_embeddings=True)
            )
        self.embeddings=output


    def align(self, other, destination):
        if self.pdb is None or other.pdb is None:
            raise ValueError("Cannot align structures without a PDB")

        command=["mustang", "-i", os.path.abspath(self.pdb), os.path.abspath(other.pdb), "-o", os.path.abspath(destination), "-r", "ON"]
        process=subprocess.run(command)
        if process.returncode != 0:
            raise ValueError("There was an error aligning structures. See error below \n {}".format(process.stderr))

        aligned_pdb=os.path.abspath(destination+".pdb")
        rotation_file=os.path.abspath(destination+".rot")
        return aligned_pdb, rotation_file


    def predict(self, sequence):
        pass


class ProteinComplex:
    def __init__(self, sequences, stoichiometry, pdb=None):
        self.sequences = sequences
        self.stoichiometry = stoichiometry
        self.pdb = pdb

    def read_pdb(self, path):
        pass

    def contacts(self):
        pass

    def predict(self, use_msa=False, model="AF3"):
        pass