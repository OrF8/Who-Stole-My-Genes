def normalize_species_name(organism: str) -> str:
    """
    Normalize organism name to "Genus species" format.
    :param organism: The input organism name.
    :return: The normalized organism name.
    """
    parts = organism.strip().split()
    return " ".join(parts[:2]) if len(parts) >= 2 else organism.strip()