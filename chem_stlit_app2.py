import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Draw
import io

# -----------------------------
# Helper Functions
# -----------------------------

def resolve_to_chembl_id(compound_name):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q={compound_name}"
    try:
        res = requests.get(url, timeout=10)
        res.raise_for_status()
        data = res.json()
        if data["page_meta"]["total_count"] > 0:
            return data["molecules"][0]["molecule_chembl_id"]
    except requests.exceptions.RequestException as e:
        st.error(f"⚠️ Failed to connect to ChEMBL for name resolution: {e}")
    return None

def get_compound_metadata(chembl_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    try:
        res = requests.get(url, timeout=10)
        res.raise_for_status()
        data = res.json()
        return {
            "Name": data.get("pref_name", "N/A"),
            "ChEMBL ID": chembl_id,
            "Molecular Weight": data.get("molecule_properties", {}).get("full_molweight", "N/A"),
            "SMILES": data.get("molecule_structures", {}).get("canonical_smiles", None),
        }
    except requests.exceptions.RequestException as e:
        st.error(f"⚠️ Failed to retrieve compound metadata: {e}")
    return None

def get_bioactivity_data(chembl_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={chembl_id}&limit=1"
    try:
        res = requests.get(url, timeout=10)
        res.raise_for_status()
        data = res.json()
        return data["page_meta"]["total_count"] > 0
    except requests.exceptions.RequestException as e:
        st.error(f"⚠️ Failed to retrieve bioactivity data: {e}")
    return False

def render_molecule_image(smiles):
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=(300, 300))
    bio = io.BytesIO()
    img.save(bio, format="PNG")
    bio.seek(0)
    return bio

# -----------------------------
# Streamlit UI
# -----------------------------

st.title("🔬 ChEMBL Compound Explorer")

compound_input = st.text_area("Enter compound names or ChEMBL IDs (one per line):")
compound_list = [c.strip() for c in compound_input.splitlines() if c.strip()]

for compound in compound_list:
    st.markdown(f"### 🔍 {compound}")
    chembl_id = resolve_to_chembl_id(compound) or (compound if compound.upper().startswith("CHEMBL") else None)

    if chembl_id:
        compound_info = get_compound_metadata(chembl_id)
        has_bioactivity = get_bioactivity_data(chembl_id)

        if compound_info:
            with st.expander(f"🧪 {compound_info['Name'] or chembl_id}"):
                # Render molecule image from SMILES using RDKit and display in Streamlit
                img_buffer = render_molecule_image(compound_info['SMILES'])
                if img_buffer:
                    st.image(img_buffer, width=300)
                else:
                    st.warning("⚠️ Unable to render molecule image.")

                st.write(f"**ChEMBL ID:** {compound_info['ChEMBL ID']}")
                st.write(f"**Molecular Weight:** {compound_info['Molecular Weight']}")
                st.write(f"**SMILES:** `{compound_info['SMILES']}`")
                st.markdown(f"[🔗 View on ChEMBL](https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/)", unsafe_allow_html=True)

                if has_bioactivity:
                    st.success("✅ Bioactivity data available.")
                else:
                    st.warning("⚠️ No bioactivity data found.")
        else:
            st.error("❌ Metadata not found.")
    else:
        st.error("❌ Compound not found in ChEMBL.")

st.markdown("[Molsoft L.L.C.: Drug-Likeness and molecular property prediction](https://molsoft.com/mprop/)")
