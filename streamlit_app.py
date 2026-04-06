import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol

st.title("🧪 Drug Molecule Analyzer (2D + 3D + Chirality)")

# --- Drug name to SMILES (basic dictionary) ---
drug_db = {
    "ketamine": "CN(C)C1=CC=CC=C1C(=O)[C@@H]2CCCC2Cl",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O"
}

# --- Input ---
drug_name = st.text_input("Enter Drug Name (optional):")
smiles_input = st.text_input("OR Enter SMILES:")

# --- Determine SMILES ---
smiles_string = None

if drug_name:
    key = drug_name.lower()
    if key in drug_db:
        smiles_string = drug_db[key]
        st.success(f"Loaded {drug_name} from database")
    else:
        st.warning("Drug not in database. Please enter SMILES.")

if smiles_input:
    smiles_string = smiles_input

# --- Process Molecule ---
if smiles_string:
    mol = Chem.MolFromSmiles(smiles_string)

    if mol is None:
        st.error("Invalid SMILES")
    else:
        # --- Chirality ---
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        st.subheader("🔬 Chiral Centers & R/S Configuration")

        if chiral_centers:
            for idx, config in chiral_centers:
                st.write(f"Atom Index: {idx} → Configuration: {config}")
        else:
            st.write("No chiral centers detected")

        # --- 2D Structure ---
        st.subheader("🖼️ 2D Structure")

        chiral_atoms = [c[0] for c in chiral_centers]

        img = Draw.MolToImage(
            mol,
            size=(400, 400),
            highlightAtoms=chiral_atoms,
            highlightColor=(1, 0, 0)
        )

        st.image(img)

        # --- 3D Structure ---
        st.subheader("🧬 3D Structure")

        mol_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv2())
        AllChem.MMFFOptimizeMolecule(mol_h)

        pdb = Chem.MolToPDBBlock(mol_h)

        viewer = py3Dmol.view(width=600, height=400)
        viewer.addModel(pdb, "pdb")
        viewer.setStyle({"stick": {}})

        for idx in chiral_atoms:
            viewer.addStyle(
                {"serial": idx + 1},
                {"sphere": {"radius": 0.6, "color": "red"}}
            )

        viewer.zoomTo()

        st.components.v1.html(viewer._make_html(), height=400)
