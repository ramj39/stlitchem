import streamlit as st
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen, QED, Lipinski, AllChem, DataStructs
import matplotlib.pyplot as plt
import subprocess
st.set_page_config(page_title="Molecular Structure & Drug-Likeness", layout="wide")
st.title("🔬 Molecular Structure & Drug-Likeness Tool")

# 🧠 Name to SMILES
def name_to_smiles(name):
    try:
        compound = pcp.get_compounds(name, 'name')
        if compound:
            return compound[0].canonical_smiles
    except:
        return None

# 🧠 SMILES to Mol
def smiles_to_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        AllChem.Compute2DCoords(mol)
    return mol

# 🧠 Molecule Image
def draw_structure(mol):
    return Draw.MolToImage(mol, size=(300, 300))

# 🧠 Property Engine with Lipinski & Toxicity
def compute_properties(mol):
    smi = Chem.MolToSmiles(mol)
    weight = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    violations = []
    if weight > 500: violations.append("MW > 500")
    if logp > 5: violations.append("LogP > 5")
    if hbd > 5: violations.append("HBD > 5")
    if hba > 10: violations.append("HBA > 10")
    lipinski = "Pass ✅" if not violations else f"Fail ❌ ({', '.join(violations)})"

    tox_flags = []
    if "N(c1ccccc1)" in smi or "[N+](=O)[O-]" in smi: tox_flags.append("Aromatic amine / Nitro")
    if any(hal in smi for hal in ["Cl", "Br", "I"]): tox_flags.append("Halogen")
    toxicity = ", ".join(tox_flags) if tox_flags else "No alerts ✅"

    return {
        "Molecular Weight": weight,
        "TPSA": Descriptors.TPSA(mol),
        "XLogP": logp,
        "H-Bond Donors": hbd,
        "H-Bond Acceptors": hba,
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "QED": QED.qed(mol),
        "Lipinski Rule of 5": lipinski,
        "Toxicity Flags": toxicity
    }

# 🧠 Similarity
def tanimoto_similarity(m1, m2):
    fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, 2)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(m2, 2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

# --- Sidebar Input ---
with st.sidebar:
    st.header("🔍 Enter Compound Names")
    user_input = st.text_area("One name per line", "Aspirin\nIbuprofen\nCaffeine")
    go_button = st.button("Analyze")
st.markdown("---")
st.markdown("Need internet-powered results?")
if st.button("🔎 Launch Similarity Search"):
    subprocess.Popen(["streamlit", "run", "compound_similarity_web.py"])
if st.sidebar.button("🧬 Launch Bioactivity Dashboard"):
    subprocess.Popen(["streamlit", "run", "compound_chembl_bioactivity.py"])
    
# 🧪 Main Logic
if go_button:
    names = [n.strip() for n in user_input.splitlines() if n.strip()]
    rows = []
    mols = []
    smiles_used = []

    for name in names:
        smi = name_to_smiles(name)
        mol = smiles_to_mol(smi) if smi else None
        if mol:
            props = compute_properties(mol)
            props["Compound"] = name
            props["SMILES"] = smi
            rows.append(props)
            mols.append(mol)
            smiles_used.append(smi)
        else:
            st.warning(f"⚠ Could not retrieve structure for: {name}")

    if rows:
        df = pd.DataFrame(rows).set_index("Compound")
        st.subheader("📋 Compound Properties (Lipinski & Toxicity)")
        st.dataframe(df)

        st.subheader("🧪 Structure Viewer")
        cols = st.columns(3)
        for i, mol in enumerate(mols):
            with cols[i % 3]:
                st.image(draw_structure(mol), caption=names[i])

        st.subheader("🧬 Tanimoto Similarity Matrix")
        sim_matrix = pd.DataFrame(index=names, columns=names)
        for i in range(len(mols)):
            for j in range(len(mols)):
                sim = tanimoto_similarity(mols[i], mols[j])
                sim_matrix.iloc[i, j] = round(sim, 3)
        st.dataframe(sim_matrix)

        st.markdown("**🔵 Similarity Heatmap**")
        fig, ax = plt.subplots()
        cax = ax.matshow(sim_matrix.astype(float), cmap="Blues")
        ax.set_xticks(range(len(names)))
        ax.set_yticks(range(len(names)))
        ax.set_xticklabels(names, rotation=90)
        ax.set_yticklabels(names)
        fig.colorbar(cax)
        st.pyplot(fig)

        # 🎯 Similarity Explorer
        st.sidebar.markdown("---")
        st.sidebar.header("🔬 Similar Molecule Finder")
        selected_name = st.sidebar.selectbox("Select reference compound", names)
        top_n = st.sidebar.slider("Top similar molecules to display", 1, len(names)-1, min(5, len(names)-1))
        threshold = st.sidebar.slider("Similarity threshold", 0.0, 1.0, 0.75, 0.01)

        st.subheader(f"🔍 Most Similar to {selected_name} (≥ {threshold:.2f})")

        idx = names.index(selected_name)
        target_mol = mols[idx]
        results = []

        for i, mol in enumerate(mols):
            if i != idx:
                score = tanimoto_similarity(target_mol, mol)
                if score >= threshold:
                    results.append((names[i], smiles_used[i], mol, score))

        results = sorted(results, key=lambda x: -x[3])[:top_n]

        if results:
            for name2, smi2, mol2, simval in results:
                col1, col2 = st.columns([1, 3])
                with col1:
                    st.image(draw_structure(mol2), width=150)
                with col2:
                    st.markdown(f"**{name2}**")
                    st.markdown(f"`{smi2}`")
                    st.markdown(f"Similarity: **{simval:.3f}**")
        else:
            st.info("No compounds matched the selected threshold.")
    else:
        st.error("❌ No valid compounds were processed.")
