import pandas as pd
from typing import Dict, Union
from ..constants import PROTEIN_ID_KEY, SEQUENCE_KEY, ORGANISM_KEY, PROTEIN_NAME_KEY, SEQ_LENGTH_KEY


def create_node_attributes(df: pd.DataFrame) -> Dict[str, Dict[str, Union[str, int]]]:
    """
    Create a dictionary of node attributes from the given DataFrame.
    :param df: The input DataFrame containing node data.
    :return: A dictionary where keys are node IDs and values are tuples of (seq_len, organism, name).
    """
    node_attributes = {}
    for _, row in df.iterrows():
        node_id = row[PROTEIN_ID_KEY]
        seq_len = len(row[SEQUENCE_KEY])
        organism = row[ORGANISM_KEY]
        name = row[PROTEIN_NAME_KEY]
        node_attributes[node_id] = {
            SEQ_LENGTH_KEY: seq_len,
            ORGANISM_KEY: organism,
            PROTEIN_NAME_KEY: name
        }
    return node_attributes


def extract_seqs(df: pd.DataFrame) -> Dict[str, str]:
    """
    Extract sequences from the DataFrame and return them as a dictionary.
    :param df: The input DataFrame containing sequence data.
    :return: A dictionary where keys are protein IDs and values are sequences.
    """
    sequences = {}
    for _, row in df.iterrows():
        protein_id = row[PROTEIN_ID_KEY]
        sequence = row[SEQUENCE_KEY]
        sequences[protein_id] = sequence.strip().upper()
    return sequences