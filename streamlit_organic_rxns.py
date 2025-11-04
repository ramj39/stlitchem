import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import MolDraw2DCairo
import io
from PIL import Image

# Set page configuration
st.set_page_config(
    page_title="Organic Chemistry Reaction Visualizer",
    page_icon="🧪",
    layout="wide"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 2rem;
    }
    .reaction-section {
        background-color: #f0f2f6;
        padding: 20px;
        border-radius: 10px;
        margin: 10px 0;
    }
    .compound-card {
        background-color: white;
        padding: 15px;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        margin: 10px 0;
    }
</style>
""", unsafe_allow_html=True)

def mol_to_image(mol, width=300, height=200):
    """Convert RDKit molecule to PIL Image"""
    if mol is None:
        return None
    drawer = MolDraw2DCairo(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png_data = drawer.GetDrawingText()
    return Image.open(io.BytesIO(png_data))

def display_compound(comp_name, smiles, description=""):
    """Display compound with name, structure, and description"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = mol_to_image(mol)
        col1, col2 = st.columns([1, 2])
        with col1:
            st.image(img, caption=comp_name, use_column_width=True)
        with col2:
            st.write(f"{comp_name}")
            st.write(f"SMILES: {smiles}")
            if description:
                st.write(description)
    else:
        st.error(f"Invalid SMILES for {comp_name}: {smiles}")

def main():
    st.markdown('<h1 class="main-header">🧪 Organic Chemistry Reaction Visualizer</h1>', 
                unsafe_allow_html=True)
    
    # Sidebar for navigation
    st.sidebar.title("Navigation")
    reaction_type = st.sidebar.selectbox(
        "Select Reaction Type",
        ["Oxidation", "Reduction", "Rearrangement", "Substitution", "Elimination", "Addition"]
    )
    
    # Compound A input
    st.sidebar.header("Compound A Input")
    custom_mode = st.sidebar.checkbox("Use custom compound")
    
    if custom_mode:
        compound_a_smiles = st.sidebar.text_input("Enter SMILES for Compound A", "CCO")
        compound_a_name = st.sidebar.text_input("Compound A Name", "Ethanol")
    else:
        predefined_compounds = {
            "Ethanol": "CCO",
            "Methanol": "CO",
            "Acetaldehyde": "CC=O",
            "Acetic Acid": "CC(=O)O",
            "Cyclohexanol": "C1CCC(CC1)O",
            "Benzaldehyde": "c1ccc(cc1)C=O",
            "2-Propanol": "CC(O)C",
            "1-Butanol": "CCCCO"
        }
        selected_compound = st.sidebar.selectbox("Select Compound A", list(predefined_compounds.keys()))
        compound_a_name = selected_compound
        compound_a_smiles = predefined_compounds[selected_compound]
    
    # Main content area
    st.header(f"{reaction_type} Reactions")
    
    # Display Compound A
    st.subheader("Starting Compound")
    display_compound(compound_a_name, compound_a_smiles)
    
    # Reaction examples based on selected type
    if reaction_type == "Oxidation":
        show_oxidation_reactions(compound_a_name, compound_a_smiles)
    elif reaction_type == "Reduction":
        show_reduction_reactions(compound_a_name, compound_a_smiles)
    elif reaction_type == "Rearrangement":
        show_rearrangement_reactions(compound_a_name, compound_a_smiles)
    elif reaction_type == "Substitution":
        show_substitution_reactions(compound_a_name, compound_a_smiles)
    elif reaction_type == "Elimination":
        show_elimination_reactions(compound_a_name, compound_a_smiles)
    elif reaction_type == "Addition":
        show_addition_reactions(compound_a_name, compound_a_smiles)

def show_oxidation_reactions(compound_name, smiles):
    """Display oxidation reactions"""
    st.markdown('<div class="reaction-section">', unsafe_allow_html=True)
    
    # Common oxidation reactions database
    oxidation_reactions = {
        "Alcohol to Aldehyde": {
            "reactant_smiles": ["CCO"],  # Ethanol
            "reactant_smiles": ["CCCO"],  # Propanol 
            "reactant_smiles": ["CCCCO"],  # Butanol
            "product_smiles": "CC=O",# Acetaldehyde
            "product_smiles": "CCC=O",# Prapanaldehyde
            "product_smiles": "CCCC=O",# Butanaldehyde
            "reagents": "PCC, CrO₃, or KMnO₄",
            "conditions": "Mild oxidation, anhydrous conditions"
        },
        "Alcohol to Carboxylic Acid": {
            "reactant_smiles": ["CCO", "CCCO"],  # Primary alcohols
            "product_smiles": "C(=O)O",  # Carboxylic acid pattern
            "reagents": "KMnO₄, K₂Cr₂O₇/H₂SO₄",
            "conditions": "Strong oxidation, acidic conditions"
        },
        "Aldehyde to Carboxylic Acid": {
            "reactant_smiles": ["CC=O"],  # Acetaldehyde
            "product_smiles": "CC(=O)O",  # Acetic acid
            "reagents": "Tollens' reagent, KMnO₄",
            "conditions": "Mild conditions"
        }
    }
    
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # Check what functional groups are present
        alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
        aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
        
        if mol.HasSubstructMatch(alcohol_pattern):
            st.success("🔍 Compound A contains alcohol functional group")
            show_reaction_example("Alcohol to Aldehyde", oxidation_reactions["Alcohol to Aldehyde"], 
                                 compound_name, smiles)
            show_reaction_example("Alcohol to Carboxylic Acid", oxidation_reactions["Alcohol to Carboxylic Acid"],
                                 compound_name, smiles)
        
        elif mol.HasSubstructMatch(aldehyde_pattern):
            st.success("🔍 Compound A contains aldehyde functional group")
            show_reaction_example("Aldehyde to Carboxylic Acid", oxidation_reactions["Aldehyde to Carboxylic Acid"],
                                 compound_name, smiles)
        else:
            st.info("💡 Try compounds like ethanol (CCO) or acetaldehyde (CC=O) for oxidation examples")
    
    st.markdown('</div>', unsafe_allow_html=True)

def show_reduction_reactions(compound_name, smiles):
    """Display reduction reactions"""
    st.markdown('<div class="reaction-section">', unsafe_allow_html=True)
    
    reduction_reactions = {
        "Aldehyde to Primary Alcohol": {
            "reactant_smiles": ["CC=O"],
            "product_smiles": "CCO",
            "reagents": "NaBH₄, LiAlH₄",
            "conditions": "Room temperature"
        },
        "Ketone to Secondary Alcohol": {
            "reactant_smiles": ["CC(=O)C"],
            "product_smiles": "CC(O)C",
            "reagents": "NaBH₄, LiAlH₄",
            "conditions": "Room temperature"
        },
        "Carboxylic Acid to Alcohol": {
            "reactant_smiles": ["CC(=O)O"],
            "product_smiles": "CCO",
            "reagents": "LiAlH₄",
            "conditions": "Anhydrous conditions"
        }
    }
    
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
        ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
        acid_pattern = Chem.MolFromSmarts("C(=O)O")
        
        if mol.HasSubstructMatch(aldehyde_pattern):
            st.success("🔍 Compound A contains aldehyde functional group")
            show_reaction_example("Aldehyde to Primary Alcohol", reduction_reactions["Aldehyde to Primary Alcohol"],
                                 compound_name, smiles)
        
        elif mol.HasSubstructMatch(ketone_pattern):
            st.success("🔍 Compound A contains ketone functional group")
            show_reaction_example("Ketone to Secondary Alcohol", reduction_reactions["Ketone to Secondary Alcohol"],
                                 compound_name, smiles)
        
        elif mol.HasSubstructMatch(acid_pattern):
            st.success("🔍 Compound A contains carboxylic acid functional group")
            show_reaction_example("Carboxylic Acid to Alcohol", reduction_reactions["Carboxylic Acid to Alcohol"],
                                 compound_name, smiles)
        else:
            st.info("💡 Try compounds like acetaldehyde (CC=O) or acetone (CC(=O)C) for reduction examples")
    
    st.markdown('</div>', unsafe_allow_html=True)

def show_rearrangement_reactions(compound_name, smiles):
    """Display rearrangement reactions"""
    st.markdown('<div class="reaction-section">', unsafe_allow_html=True)
    
    rearrangement_examples = {
        "Pinacol Rearrangement": {
            "reactant_smiles": "CC(C)(C)C(C)(C)O",
            "product_smiles": "CC(C)(C)C(=O)C",
            "reagents": "H₂SO₄, H₃O⁺",
            "conditions": "Acid-catalyzed",
            "description": "Rearrangement of 1,2-diols to carbonyl compounds"
        },
        "Beckmann Rearrangement": {
            "reactant_smiles": "CC(=NOH)C",
            "product_smiles": "CC(=O)NC",
            "reagents": "Acid catalyst",
            "conditions": "Acidic conditions",
            "description": "Rearrangement of oximes to amides"
        },
        "Hofmann Rearrangement": {
            "reactant_smiles": "CC(=O)N",
            "product_smiles": "CNC=O",
            "reagents": "Br₂, NaOH",
            "conditions": "Basic conditions",
            "description": "Conversion of primary amides to amines with one less carbon"
        }
    }
    
    st.subheader("Common Rearrangement Reactions")
    
    for reaction_name, data in rearrangement_examples.items():
        with st.expander(f"🎯 {reaction_name}"):
            st.write(f"*Description:* {data['description']}")
            st.write(f"*Reagents:* {data['reagents']}")
            st.write(f"*Conditions:* {data['conditions']}")
            
            col1, col2, col3 = st.columns([1, 1, 1])
            with col1:
                display_compound("Reactant", data["reactant_smiles"])
            with col2:
                st.markdown("<h3 style='text-align: center;'>→</h3>", unsafe_allow_html=True)
                st.write("*Rearrangement*")
            with col3:
                display_compound("Product", data["product_smiles"])
    
    st.markdown('</div>', unsafe_allow_html=True)

def show_substitution_reactions(compound_name, smiles):
    """Display substitution reactions"""
    st.markdown('<div class="reaction-section">', unsafe_allow_html=True)
    
    substitution_examples = {
        "SN2 Reaction": {
            "reactant_smiles": "CCl",
            "product_smiles": "CI",
            "reagents": "NaI in acetone",
            "conditions": "Polar aprotic solvent",
            "description": "Bimolecular nucleophilic substitution"
        },
        "SN1 Reaction": {
            "reactant_smiles": "C(C)(C)Cl",
            "product_smiles": "C(C)(C)O",
            "reagents": "H₂O",
            "conditions": "Polar protic solvent",
            "description": "Unimolecular nucleophilic substitution"
        }
    }
    
    st.subheader("Substitution Reactions")
    st.info("SN1 and SN2 mechanisms")
    
    for reaction_name, data in substitution_examples.items():
        with st.expander(f"🎯 {reaction_name}"):
            st.write(f"*Description:* {data['description']}")
            st.write(f"*Reagents:* {data['reagents']}")
            st.write(f"*Conditions:* {data['conditions']}")
            
            col1, col2, col3 = st.columns([1, 1, 1])
            with col1:
                display_compound("Reactant", data["reactant_smiles"])
            with col2:
                st.markdown("<h3 style='text-align: center;'>→</h3>", unsafe_allow_html=True)
                st.write("*Substitution*")
            with col3:
                display_compound("Product", data["product_smiles"])
    
    st.markdown('</div>', unsafe_allow_html=True)

def show_elimination_reactions(compound_name, smiles):
    """Display elimination reactions"""
    st.markdown('<div class="reaction-section">', unsafe_allow_html=True)
    
    elimination_examples = {
        "E2 Elimination": {
            "reactant_smiles": "CCCCl",
            "product_smiles": "C=C",
            "reagents": "KOH, ethanol",
            "conditions": "Strong base, heat",
            "description": "Bimolecular elimination"
        },
        "E1 Elimination": {
            "reactant_smiles": "CC(C)(C)Cl",
            "product_smiles": "C=C(C)C",
            "reagents": "H₂O, heat",
            "conditions": "Weak base, thermal",
            "description": "Unimolecular elimination"
        }
    }
    
    st.subheader("Elimination Reactions")
    st.info("E1 and E2 mechanisms forming alkenes")
    
    for reaction_name, data in elimination_examples.items():
        with st.expander(f"🎯 {reaction_name}"):
            st.write(f"*Description:* {data['description']}")
            st.write(f"*Reagents:* {data['reagents']}")
            st.write(f"*Conditions:* {data['conditions']}")
            
            col1, col2, col3 = st.columns([1, 1, 1])
            with col1:
                display_compound("Reactant", data["reactant_smiles"])
            with col2:
                st.markdown("<h3 style='text-align: center;'>→</h3>", unsafe_allow_html=True)
                st.write("*Elimination*")
            with col3:
                display_compound("Product", data["product_smiles"])
    
    st.markdown('</div>', unsafe_allow_html=True)

def show_addition_reactions(compound_name, smiles):
    """Display addition reactions"""
    st.markdown('<div class="reaction-section">', unsafe_allow_html=True)
    
    addition_examples = {
        "Hydrogenation": {
            "reactant_smiles": "C=C",
            "product_smiles": "CC",
            "reagents": "H₂, Pt/Pd/Ni catalyst",
            "conditions": "Room temperature, pressure",
            "description": "Addition of hydrogen to alkenes/alkynes"
        },
        "Halogen Addition": {
            "reactant_smiles": "C=C",
            "product_smiles": "CC(Br)Br",
            "reagents": "Br₂ in CCl₄",
            "conditions": "Room temperature",
            "description": "Anti addition of halogens"
        }
    }
    
    st.subheader("Addition Reactions")
    st.info("Electrophilic and nucleophilic addition to multiple bonds")
    
    for reaction_name, data in addition_examples.items():
        with st.expander(f"🎯 {reaction_name}"):
            st.write(f"*Description:* {data['description']}")
            st.write(f"*Reagents:* {data['reagents']}")
            st.write(f"*Conditions:* {data['conditions']}")
            
            col1, col2, col3 = st.columns([1, 1, 1])
            with col1:
                display_compound("Reactant", data["reactant_smiles"])
            with col2:
                st.markdown("<h3 style='text-align: center;'>→</h3>", unsafe_allow_html=True)
                st.write("*Addition*")
            with col3:
                display_compound("Product", data["product_smiles"])
    
    st.markdown('</div>', unsafe_allow_html=True)

def show_reaction_example(reaction_name, reaction_data, compound_name, compound_smiles):
    """Display a specific reaction example"""
    with st.expander(f"🎯 {reaction_name}"):
        st.write(f"*Reagents:* {reaction_data['reagents']}")
        st.write(f"*Conditions:* {reaction_data['conditions']}")
        
        col1, col2, col3 = st.columns([1, 1, 1])
        
        with col1:
            display_compound("Reactant", compound_smiles, compound_name)
        
        with col2:
            st.markdown("<h3 style='text-align: center;'>→</h3>", unsafe_allow_html=True)
            st.write("*Reaction*")
        
        with col3:
            # Try to generate product based on reaction pattern
            product_smiles = predict_product(compound_smiles, reaction_name)
            if product_smiles:
                display_compound("Product", product_smiles, "Predicted Product")
            else:
                # Show example product if prediction not available
                if "product_smiles" in reaction_data:
                    display_compound("Example Product", reaction_data["product_smiles"], "Example Reaction")
                else:
                    st.info("Product prediction not available for this compound")

def predict_product(smiles, reaction_type):
    """Simple product prediction based on reaction type"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        if "Alcohol to Aldehyde" in reaction_type:
            # Simple transformation: primary alcohol to aldehyde
            pattern = Chem.MolFromSmarts("[CH2][OH]")
            if mol.HasSubstructMatch(pattern):
                # This is a simplified transformation
                return smiles.replace("CO", "C=O")
        
        elif "Alcohol to Carboxylic Acid" in reaction_type:
            pattern = Chem.MolFromSmarts("[CH2][OH]")
            if mol.HasSubstructMatch(pattern):
                return smiles.replace("CO", "C(=O)O")
        
        elif "Aldehyde to Carboxylic Acid" in reaction_type:
            pattern = Chem.MolFromSmarts("C=O")
            if mol.HasSubstructMatch(pattern):
                return smiles.replace("C=O", "C(=O)O")
    
    return None

# Add information section
def add_info_section():
    st.sidebar.markdown("---")
    st.sidebar.header("About")
    st.sidebar.info("""
    This application demonstrates organic chemistry reactions using RDKit for structure visualization.
    
    *Features:*
    - View compound structures
    - Explore different reaction types
    - See reaction conditions and reagents
    - Custom compound input via SMILES
    """)

if __name__ == "__main__":
    add_info_section()
    main()
