

class MoleculeBinder:
    def __init__(self, protein, pockets=None, cleanup=True):
        """
        the structure needed for the molecule binder
        :param protein: a Protein object, we need a structure, if the protein does not have a pdb file we can
        either use the structure module to predict one or download from pdb doesn't matter
        :param pockets: list of pocket coordinates, if none we can use find pockets method
        :param cleanup: whether to cleanup the pdb, that is remove non protein molecules waters etc. This will use pdbcleaner
        """
        pass

    def find_pockets(self, cleanup=True, return_pockets=True, destination=None):
        """
        this runs fpocket and gets the pocket pdb as well as pocket properties in a dict
        :param cleanup: cleanup the intermediate files.
        :param return_pockets: whether to return pockets or not as a dict if false will just keep the files
        :param destination: where to save potentially intermediate files
        :return: message and dict
        """
        pass

    def define_pockets(self, aas):
        """
        for a given set of aas define a pocket
        :param aas: set of aas to consider
        :return: a list of coordinates that capture the aas specified
        """
        pass
    def generate_molecules(self, pocket, n):
        """
        this runs pocket2mol for a give pocket and generated n molecules. The pockets can be de
        :param pocket:
        :param n:
        :return:
        """
        pass

    def search_molecules(self, library, using="ecfp4", metric="tanimoto", n=100):
        """
        using useard molecules will search library for a given metric, the only metric implemented is tanimoto similarity
        :param library: a list of directories that contain usearch parquet files
        :param using: the implemented ones are ecfp4, fcfp4 and maccs
        :param metric: tanimoto otherwise NotImplementedError
        :param n: number of molecules to return
        :return: list of smiles for molecules
        """
        pass

    def calculate_properties(self, *args):
        """
        calculate numerous properties using RDkit simple wrapper
        :param args: properties to get
        :return: returns properites of a smiles string based on args
        """

    def verify(self, method="AF3", **kwargs):
        """
        run AF3 to re-predict the structure to verify binding
        :param method:
        :param kwargs:
        :return:
        """
        pass


class PeptideBinder:
    def __init__(self, protein, region):
        self.protein = protein
        if region is not None:
            self.region = region

    def find_binding_site(self):
        """
        this will again use fpocket
        :return:
        """
        pass

    def generate_peptides(self, region, n):
        pass

    def verify(self, method="AF3", **kwargs):
        pass