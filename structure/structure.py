import biotite
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig

from ccm_demo.simulation import simulation

class Structure:
    def __init__(self, pdb):
        self.pdb = pdb
        self.sasa=None,
        self.embeddings=None
        self.foldseek=None
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to(self.device)

    def download(self):
        pass

    def simulate(self, simulation):
        pass

    def embeddings(self, model="esm3", normalize=False):
        pass

    def align(self, other):
        pass

    def predict(self, sequence):
        pass


class ProteinComplex:
    def __init__(self, pdb):
        self.pdb = pdb

    def contacts(self):
        pass

    def predict(self, sequences, stoichiometry):
        pass