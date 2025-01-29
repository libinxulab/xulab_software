"""
    toxcccs/utils/smiles.py
    
    Ryan Nguyen (ryan97@uw.edu)
    12/20/24

    description:
            Module with utiliy for converting SMILES strings to chemical structures.
"""

from rdkit import Chem
from rdkit.Chem import Draw


def smiles_to_structure(smiles, img_size=(100, 100)):
    """
    utils/smiles
    description:
            Converts a SMILES string to a chemical structure.
    parameters:
            smiles (str) -- SMILES string.
            img_ize (tuple) -- image size. Default is (100, 100).
    returns:
            (PIL.Image.Image or None) -- Image of the chemical structure.
    """

    # Check if the string input is valid
    if not isinstance(smiles, str) or smiles.strip() == "":
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=img_size)
        return img
    else:
        return None
