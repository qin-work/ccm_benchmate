import json
import string

from ccm_demo.containers.utils import *


def generate_alpha_id(n):
    alpha_IDs = list(string.ascii_uppercase)
    if n > 26:  # Check if AA-ZZ should be appended to alpha_IDs
        for first in string.ascii_uppercase:
            for second in string.ascii_uppercase:
                alpha_IDs.append(first + second)
                if len(alpha_IDs) >= n:  # Ensures unnecessary IDs are not generated
                    return alpha_IDs
    return alpha_IDs[:n]

#TODO needs to account for small molecules, DNA, RNA
def generate_json(name, sequences, stoichiometry, fpath, seed=1, version=1):
    total_chains = sum(map(int, stoichiometry))
    alpha_IDs = generate_alpha_id(total_chains)

    dialect = "alphafold3"
    version = version
    if len(sequences) != len(stoichiometry):
        raise ValueError("sequences and stoichiometry must have same length")

    protein_sequences = []
    alpha_index = 0

    for seq, stoich in zip(sequences, stoichiometry):
        ids = alpha_IDs[alpha_index:alpha_index + stoich]
        id_value = ids[0] if len(ids) == 1 else ids  # Convert to string if one item
        protein_sequences.append({"protein": {"id": id_value, "sequence": seq}})
        alpha_index += stoich

    json_file = {"name": name, "modelSeeds": [412], "sequences": protein_sequences, "dialect": dialect, "version": 1}
    file_name = f"{name}.json"
    json_path = os.path.abspath(os.path.join(fpath, file_name))
    with open(json_path, "w") as f:
        json.dump(json_file, f)

# this is for demo purposes
command_dict = {"bind_mounts": [{"local": "/hpf/largeprojects/ccmbio/acelik_files/alphafold3",
                                  "target": "/af3"},
                                 {"local":None,
                                  "target":"/inputs"},
                                 {"local":None,
                                  "target":"/outputs"}],
                "singularity_args":["--nv", "--cleanenv"],
                 "run_command":"""
                 python /af3/run_alphafold.py \\
                    --db_dir /af3/af3_db \\
                    --input_dir /inputs \\
                    --json_path /inputs/$file \\
                    --model_dir /af3/af3_weights \\
                    --output_dir /outputs \\
                    --run_data_pipeline={} \\
                    --run_inference={}
                 """}
