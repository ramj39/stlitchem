import streamlit as st
from chempy import balance_stoichiometry
from chempy.util import periodic
import re

# Initialize session state for history
if 'history' not in st.session_state:
    st.session_state.history = []

# Expanded chemical dictionary
COMMON_NAMES = {
    "benzaldehyde": "C6H5CHO", "benzene": "C6H6", "toluene": "C6H5CH3", 
    "phenol": "C6H5OH", "aniline": "C6H5NH2", "benzoic acid": "C6H5COOH",
    "methanol": "CH3OH", "ethanol": "C2H5OH", "formic acid": "HCOOH", 
    "acetic acid": "CH3COOH", "oxygen": "O2", "hydrogen peroxide": "H2O2",
    "potassium permanganate": "KMnO4", "tollens reagent": "Ag(NH3)2OH",
    "water": "H2O", "sodium hydroxide": "NaOH", "sulfuric acid": "H2SO4",
    "nitric acid": "HNO3", "hydrochloric acid": "HCl", "carbon dioxide": "CO2",
    "formaldehyde": "HCHO", "acetaldehyde": "CH3CHO", "glucose": "C6H12O6",
    "acetone": "CH3COCH3", "ethylene": "C2H4", "acetylene": "C2H2",
    "ammonia": "NH3", "bromine": "Br2", "chlorine": "Cl2"
}

# Expanded reaction database
COMMON_REACTIONS = {
    # Oxidation Reactions
    "oxidation of benzaldehyde": ("C6H5CHO + [O]", "C6H5COOH"),
    "oxidation of formaldehyde": ("HCHO + [O]", "HCOOH"),
    "oxidation of acetaldehyde": ("CH3CHO + [O]", "CH3COOH"),
    "oxidation with kmno4": ("RCHO + KMnO4 + H2SO4", "RCOOH + MnSO4 + K2SO4 + H2O"),
    "tollens test": ("RCHO + Ag(NH3)2OH", "RCOONH4 + Ag + H2O + NH3"),
    "fehling test": ("RCHO + Cu2+", "RCOOH + Cu2O"),
    
    # Named Organic Synthesis Reactions
    "cannizzaro reaction": ("2C6H5CHO", "C6H5COOH + C6H5CH2OH"),
    "aldol condensation": ("2CH3CHO", "CH3CH(OH)CH2CHO"),
    "benzoin condensation": ("2C6H5CHO", "C6H5CH(OH)C(O)C6H5"),
    "grignard reaction": ("RMgBr + R'CHO", "R-CH(OH)-R'"),
    "wittig reaction": ("R2C=O + Ph3P=CR2", "R2C=CR2 + Ph3PO"),
    "diels-alder reaction": ("diene + dienophile", "cyclohexene"),
    "fischer esterification": ("RCOOH + R'OH", "RCOOR' + H2O"),
    "claisen condensation": ("2RCH2COOR'", "RCH2C(O)CHRCOR' + R'OH"),
    "friedel-crafts acylation": ("C6H6 + RCOCl", "C6H5COR + HCl"),
    "wurtz reaction": ("2R-X + 2Na", "R-R + 2NaX"),
    "reimer-tiemann reaction": ("C6H5OH + CHCl3", "o-HOC6H4CHO"),
    "clemmensen reduction": ("R-C=O", "R-CH2"),
    "wolff-kishner reduction": ("R-C=O", "R-CH2"),
    "darzens glycidic ester condensation":("R-CO","RR-COR"),
    
    # Previous reactions
    "chlorination of benzene": ("C6H6 + Cl2", "C6H5Cl + HCl"),
    "bromination of benzene": ("C6H6 + Br2", "C6H5Br + HBr"),
    "nitration of benzene": ("C6H6 + HNO3", "C6H5NO2 + H2O"),
    "sulfonation of benzene": ("C6H6 + H2SO4", "C6H5SO3H + H2O"),
    "combustion of benzene": ("C6H6 + O2", "CO2 + H2O"),
    "combustion of ethanol": ("C2H5OH + O2", "CO2 + H2O"),
    "neutralization": ("HCl + NaOH", "NaCl + H2O"),
    "ammonia synthesis": ("N2 + H2", "NH3"),
    "photosynthesis": ("CO2 + H2O", "C6H12O6 + O2"),
    "fermentation": ("C6H12O6", "C2H5OH + CO2"),
}

# Organic synthesis reactions with descriptions
ORGANIC_SYNTHESIS = {
    "Grignard Reaction": {
        "reactants": "RMgBr + R'CHO",
        "products": "R-CH(OH)-R'",
        "description": "Nucleophilic addition of organomagnesium compounds to carbonyl groups to form alcohols",
        "mechanism": "Nucleophilic addition",
        "year": 1900,
        "chemist": "Victor Grignard"
    },
    "Wittig Reaction": {
        "reactants": "R2C=O + Ph3P=CR2",
        "products": "R2C=CR2 + Ph3PO",
        "description": "Conversion of carbonyls to alkenes using phosphonium ylides",
        "mechanism": "Nucleophilic addition then elimination",
        "year": 1954,
        "chemist": "Georg Wittig"
    },
    "Diels-Alder Reaction": {
        "reactants": "diene + dienophile",
        "products": "cyclohexene derivative",
        "description": "Cycloaddition between a conjugated diene and a substituted alkene",
        "mechanism": "Pericyclic [4+2] cycloaddition",
        "year": 1928,
        "chemist": "Otto Diels & Kurt Alder"
    },
    "Friedel-Crafts Acylation": {
        "reactants": "C6H6 + RCOCl",
        "products": "C6H5COR + HCl",
        "description": "Electrophilic aromatic substitution using acyl chlorides to form ketones",
        "mechanism": "Electrophilic aromatic substitution",
        "year": 1877,
        "chemist": "Charles Friedel & James Crafts"
    },
    "Claisen Condensation": {
        "reactants": "2RCH2COOR'",
        "products": "RCH2C(O)CHRCOR' + R'OH",
        "description": "Base-catalyzed condensation of esters to form β-ketoesters",
        "mechanism": "Nucleophilic acyl substitution",
        "year": 1887,
        "chemist": "Rainer Ludwig Claisen"
    },
    "Wurtz Reaction": {
        "reactants": "2R-X + 2Na",
        "products": "R-R + 2NaX",
        "description": "Coupling reaction of alkyl halides to form symmetric alkanes",
        "mechanism": "Radical mechanism",
        "year": 1855,
        "chemist": "Charles Adolphe Wurtz"
    },
    "Reimer-Tiemann Reaction": {
        "reactants": "C6H5OH + CHCl3",
        "products": "o-HOC6H4CHO",
        "description": "Ortho-formylation of phenols under basic conditions",
        "mechanism": "Carbene intermediate",
        "year": 1876,
        "chemist": "Karl Reimer & Ferdinand Tiemann"
    },
    "Cannizzaro Reaction": {
        "reactants": "2ArCHO",
        "products": "ArCOOH + ArCH2OH",
        "description": "Disproportionation of aldehydes lacking alpha hydrogens",
        "mechanism": "Hydride transfer",
        "year": 1853,
        "chemist": "Stanislao Cannizzaro"
    },
    "Aldol Condensation": {
        "reactants": "2RCH2CHO",
        "products": "RCH2CH(OH)CH(R)CHO",
        "description": "Formation of β-hydroxy carbonyl compounds from aldehydes/ketones",
        "mechanism": "Enolate formation then nucleophilic addition",
        "year": 1872,
        "chemist": "Charles Adolphe Wurtz"
    }
}

# Resource links for each reaction type
REACTION_RESOURCES = {
    "Benzaldehyde Oxidation": "https://en.wikipedia.org/wiki/Benzaldehyde#Reactions",
    "Cannizzaro Reaction": "https://en.wikipedia.org/wiki/Cannizzaro_reaction",
    "Aldol Condensation": "https://en.wikipedia.org/wiki/Aldol_condensation",
    "Oxidation with KMnO₄": "https://en.wikipedia.org/wiki/Potassium_permanganate#Organic_chemistry",
    "Tollens Test": "https://en.wikipedia.org/wiki/Tollens%27_reagent",
    "Fehling Test": "https://en.wikipedia.org/wiki/Fehling%27s_solution",
    "Aromatic Halogenation": "https://en.wikipedia.org/wiki/Electrophilic_halogenation",
    "Aromatic Nitration": "https://en.wikipedia.org/wiki/Nitration",
    "Aromatic Sulfonation": "https://en.wikipedia.org/wiki/Sulfonation",
    "Combustion": "https://en.wikipedia.org/wiki/Combustion",
    "Organic Reaction": "https://en.wikipedia.org/wiki/Organic_reaction",
    "Acid-Base Neutralization": "https://en.wikipedia.org/wiki/Neutralization_(chemistry)",
    "Ammonia Synthesis": "https://en.wikipedia.org/wiki/Haber_process",
    "Photosynthesis": "https://en.wikipedia.org/wiki/Photosynthesis",
    "Fermentation": "https://en.wikipedia.org/wiki/Fermentation",
    "Grignard Reaction": "https://en.wikipedia.org/wiki/Grignard_reaction",
    "Wittig Reaction": "https://en.wikipedia.org/wiki/Wittig_reaction",
    "Diels-Alder Reaction": "https://en.wikipedia.org/wiki/Diels%E2%80%93Alder_reaction",
    "Friedel-Crafts Acylation": "https://en.wikipedia.org/wiki/Friedel%E2%80%93Crafts_reaction",
    "Claisen Condensation": "https://en.wikipedia.org/wiki/Claisen_condensation",
    "Wurtz Reaction": "https://en.wikipedia.org/wiki/Wurtz_reaction",
    "Reimer-Tiemann Reaction": "https://en.wikipedia.org/wiki/Reimer%E2%80%93Tiemann_reaction"
}

# Improved formula parser
def parse_formula(formula_str):
    """Handle complex organic formulas with special groups"""
    # Replace common organic groups
    replacements = {
        'R-': 'H-', 
        'R': 'H',
        'Ph': 'C6H5',
        'Bz': 'C6H5CO',
        'Ac': 'CH3CO',
        '[O]': 'O',
        ' ': ''
    }
    
    for key, value in replacements.items():
        formula_str = formula_str.replace(key, value)
    
    return formula_str

# Enhanced chemical parser with coefficient handling
def parse_chemicals(input_str):
    input_str = input_str.lower()
    
    # Check for exact reaction matches first
    for name in COMMON_REACTIONS:
        if name in input_str:
            return COMMON_REACTIONS[name]
    
    # Tokenize and identify chemicals
    tokens = re.findall(r'\b[\w\s-]+\b', input_str)
    chemicals = []
    
    for token in tokens:
        token = token.strip()
        if token in COMMON_NAMES:
            formula = parse_formula(COMMON_NAMES[token])
            chemicals.append(formula)
        else:
            # Handle formulas with coefficients
            if re.match(r'^\d+[A-Z]', token):
                # Extract coefficient and formula
                coeff = re.match(r'^(\d+)(.*)', token).group(1)
                formula = re.match(r'^(\d+)(.*)', token).group(2)
                # Add multiple instances instead of coefficients
                formula = parse_formula(formula)
                chemicals.extend([formula] * int(coeff))
            elif any(c.isalpha() for c in token):
                # Handle complex formulas
                try:
                    parsed = parse_formula(token)
                    if all(part.capitalize() in periodic.symbols 
                           for part in re.findall('[A-Z][^A-Z]*', parsed) 
                           if part):
                        chemicals.append(parsed)
                except:
                    continue
                    
    return chemicals

# Enhanced balancing function with oxidation handling
def balance_equation(reactants_str, products_str):
    try:
        # Split and parse each component
        reactants = []
        for r in reactants_str.split('+'):
            r = r.strip()
            if re.match(r'^\d+[A-Z]', r):
                coeff = re.match(r'^(\d+)(.*)', r).group(1)
                formula = parse_formula(re.match(r'^(\d+)(.*)', r).group(2))
                reactants.extend([formula] * int(coeff))
            else:
                reactants.append(parse_formula(r))
                
        products = []
        for p in products_str.split('+'):
            p = p.strip()
            if re.match(r'^\d+[A-Z]', p):
                coeff = re.match(r'^(\d+)(.*)', p).group(1)
                formula = parse_formula(re.match(r'^(\d+)(.*)', p).group(2))
                products.extend([formula] * int(coeff))
            else:
                products.append(parse_formula(p))
        
        # Balance stoichiometry
        reac, prod = balance_stoichiometry(reactants, products)
        
        # Format equation
        reactant_side = " + ".join([f"{coeff} {formula}" if coeff != 1 else formula 
                                   for formula, coeff in reac.items()])
        product_side = " + ".join([f"{coeff} {formula}" if coeff != 1 else formula 
                                  for formula, coeff in prod.items()])
        
        # Replace H with R for generic formulas
        final_eq = f"{reactant_side} → {product_side}"
        final_eq = final_eq.replace("H-", "R-").replace("HCHO", "RCHO").replace("HCOOH", "RCOOH")
        return final_eq
    
    except Exception as e:
        return f"Error: {str(e)}"

# Enhanced reaction type detection
def detect_reaction_type(equation):
    equation = equation.lower()
    if "kmno4" in equation or "mnso4" in equation:
        return "Oxidation with KMnO₄"
    if "ag" in equation and "nh3" in equation:
        return "Tollens Test"
    if "cu2+" in equation or "cu2o" in equation:
        return "Fehling Test"
    if "c6h5cho" in equation and "c6h5cooh" in equation:
        return "Benzaldehyde Oxidation"
    if "c6h5cho" in equation and "c6h5ch2oh" in equation:
        return "Cannizzaro Reaction"
    if "ch3cho" in equation and "ch3chohch2cho" in equation:
        return "Aldol Condensation"
    if "c6h6" in equation and any(h in equation for h in ["cl2", "br2"]):
        return "Aromatic Halogenation"
    if "c6h6" in equation and "hno3" in equation:
        return "Aromatic Nitration"
    if "c6h6" in equation and "h2so4" in equation:
        return "Aromatic Sulfonation"
    if "co2" in equation and "h2o" in equation and "o2" in equation:
        return "Combustion"
    if "c6h12o6" in equation and "c2h5oh" in equation:
        return "Fermentation"
    if "h2o" in equation and "acid" in equation:
        return "Acid-Base Neutralization"
    if "nh3" in equation:
        return "Ammonia Synthesis"
    if "c6h12o6" in equation and "o2" in equation:
        return "Photosynthesis"
    if "ph3p" in equation:
        return "Wittig Reaction"
    if "rmgbr" in equation:
        return "Grignard Reaction"
    if "diene" in equation and "dienophile" in equation:
        return "Diels-Alder Reaction"
    if "rcoor" in equation:
        return "Fischer Esterification"
    return "Organic Reaction"

# Streamlit UI
st.title("🧪 Organic Synthesis Explorer")
st.subheader("Named reactions, synthesis planning, and equation balancing")

# Create tabs for different sections
tab1, tab2, tab3 = st.tabs(["Reaction Balancer", "Organic Synthesis", "Resources"])

with tab1:
    # Input section
    query = st.text_input("Ask a question about chemical reactions:", 
                          placeholder="e.g., oxidation of benzaldehyde, Cannizzaro reaction")

    # Create list of reaction types
    reaction_types = sorted({name.split()[0] for name in COMMON_REACTIONS})
    reaction_type = st.selectbox("Or choose a reaction type:", ["All"] + reaction_types)

    # Display common reactions if type selected
    if reaction_type != "All":
        st.write(f"*Common {reaction_type} reactions:*")
        for name, reaction_data in COMMON_REACTIONS.items():
            if name.startswith(reaction_type):
                if len(reaction_data) == 2:
                    react, prod = reaction_data
                    balanced_eq = balance_equation(react, prod)
                    st.code(f"{name.capitalize()}: {balanced_eq}")

    # Process question
    if st.button("Get Reaction") and query:
        with st.spinner("Balancing equation..."):
            # Try to match common reactions first
            query_lower = query.lower()
            matched = False
            
            for name, reaction_data in COMMON_REACTIONS.items():
                if name in query_lower and len(reaction_data) == 2:
                    react, prod = reaction_data
                    balanced_eq = balance_equation(react, prod)
                    rtype = detect_reaction_type(balanced_eq)
                    
                    st.success(f"{name.capitalize()} Reaction**")
                    st.subheader("Balanced Equation:")
                    st.code(balanced_eq)
                    st.write(f"*Type:* {rtype}")
                    
                    # Add resource link
                    if rtype in REACTION_RESOURCES:
                        st.markdown(f"*Learn more:* [Wikipedia article on {rtype}]({REACTION_RESOURCES[rtype]})")
                    
                    # Add to history
                    st.session_state.history.insert(0, {
                        "query": query,
                        "equation": balanced_eq,
                        "type": rtype
                    })
                    matched = True
                    break
            
            # If no common match found, try to parse
            if not matched:
                chemicals = parse_chemicals(query)
                if isinstance(chemicals, tuple) and len(chemicals) == 2:
                    react, prod = chemicals
                    balanced_eq = balance_equation(react, prod)
                    rtype = detect_reaction_type(balanced_eq)
                    
                    st.success("Balanced Equation")
                    st.subheader("Balanced Equation:")
                    st.code(balanced_eq)
                    st.write(f"*Type:* {rtype}")
                    
                    # Add resource link
                    if rtype in REACTION_RESOURCES:
                        st.markdown(f"*Learn more:* [Wikipedia article on {rtype}]({REACTION_RESOURCES[rtype]})")
                    
                    # Add to history
                    st.session_state.history.insert(0, {
                        "query": query,
                        "equation": balanced_eq,
                        "type": rtype
                    })
                elif isinstance(chemicals, list) and len(chemicals) >= 2:
                    # Try to split into reactants/products
                    reactants = " + ".join(chemicals[:len(chemicals)//2])
                    products = " + ".join(chemicals[len(chemicals)//2:])
                    
                    balanced_eq = balance_equation(reactants, products)
                    rtype = detect_reaction_type(balanced_eq)
                    
                    st.success("Balanced Equation")
                    st.subheader("Balanced Equation:")
                    st.code(balanced_eq)
                    st.write(f"*Type:* {rtype}")
                    
                    # Add resource link
                    if rtype in REACTION_RESOURCES:
                        st.markdown(f"*Learn more:* [Wikipedia article on {rtype}]({REACTION_RESOURCES[rtype]})")
                    
                    # Add to history
                    st.session_state.history.insert(0, {
                        "query": query,
                        "equation": balanced_eq,
                        "type": rtype
                    })
                else:
                    st.warning("Couldn't identify the reaction. Try examples like:")
                    st.write("- Oxidation of benzaldehyde")
                    st.write("- Cannizzaro reaction")
                    st.write("- Aldol condensation")
                    st.write("- Tollens test")

    # Display history
    # Display history - FIXED VERSION
if st.session_state.history:
    st.divider()
    st.subheader("Recent Queries")
    for i, history_item in enumerate(st.session_state.history[:5]):
        st.write(f"{i+1}. *{history_item['query']}*")
        st.caption(f"Equation: {history_item['equation']}")
        st.caption(f"Type: {history_item['type']}")
        
        if history_item['type'] in REACTION_RESOURCES:
            url = REACTION_RESOURCES[history_item['type']]
            st.caption(f"Resource: [Learn more about {history_item['type']}]({url})")
        
        st.write("")
    #if st.session_state.history:
   #     st.divider()
     #   st.subheader("Recent Queries")
       # for i, item in enumerate(st.session_state.history[:5]):  # Fixed this line
         #   st.write(f"{i+1}. *{item['query']}*")
           # st.caption(f"Equation: {item['equation']}")
            #st.caption(f"Type: {item['type']}")
            #if item['type'] in REACTION_RESOURCES:
           #     st.caption(f"Resource: [Learn more about {item['type']}]({REACTION_RESOURCES[item['type']})")
   # if item['type'] in REACTION_RESOURCES:
     #   url = REACTION_RESOURCES[item['type']]
     #   st.caption(f"Resource: [Learn more about {item['type']}]({url})")    
    #st.write("")  # This line is now fixed

with tab2:
    st.header("Organic Synthesis Reactions")
    st.subheader("Explore well-known named reactions in organic chemistry")
    
    # Search and filter
    search_term = st.text_input("Search for a named reaction:", placeholder="e.g., Grignard, Diels-Alder")
    st.caption("Browse through the most important organic synthesis reactions")
    
    # Filter reactions based on search term
    filtered_reactions = ORGANIC_SYNTHESIS
    if search_term:
        search_term = search_term.lower()
        filtered_reactions = {
            name: data for name, data in ORGANIC_SYNTHESIS.items()
            if search_term in name.lower() or 
            any(search_term in str(val).lower() for val in data.values())
        }
    
    # Display reactions in a grid
    cols = st.columns(3)
    for i, (name, data) in enumerate(filtered_reactions.items()):
        with cols[i % 3]:
            with st.expander(f"{name}"):
                st.write(f"*Reaction:* {data['reactants']} → {data['products']}")
                st.write(f"*Description:* {data['description']}")
                st.write(f"*Mechanism:* {data['mechanism']}")
                st.write(f"*Discovered by:* {data['chemist']} ({data['year']})")
                
                # Try to balance the equation
                try:
                    balanced_eq = balance_equation(data['reactants'], data['products'])
                    st.code(f"Balanced: {balanced_eq}")
                except:
                    st.warning("Couldn't balance this generic reaction")
                
                # Add resource link
                if name in REACTION_RESOURCES:
                    st.markdown(f"[Learn more about {name}]({REACTION_RESOURCES[name]})")

with tab3:
    st.header("Educational Resources")
    st.subheader("Comprehensive learning materials for organic chemistry")
    
    st.markdown("### Reaction Mechanisms")
    st.markdown("- [Organic Chemistry Portal](https://www.organic-chemistry.org/)")
    st.markdown("- [Master Organic Chemistry](https://www.masterorganicchemistry.com/)")
    st.markdown("- [Khan Academy Organic Chemistry](https://www.khanacademy.org/science/organic-chemistry)")
    st.markdown("- [Organic Chemistry Portal](https://archive.org/details/finarorganicchemistryvol1/)")
    st.markdown("### Synthesis Planning")
    st.markdown("- [SynArchive Database](https://synarchive.com/)")
    st.markdown("- [Organic Synthesis Search](https://www.organic-synthesis.org/)")
    st.markdown("- [Retrosynthesis Tools](https://scifinder-n.cas.org/)")
    st.markdown("- [organic synthesis](https://en.wikipedia.org/wiki/Organic_synthesis)")
    st.markdown("-------------------------------------------------")
    st.markdown("below are listed some named synthesis/reactions")
    st.markdown("- [Darzwns Reaction](https://en.wikipedia.org/wiki/Darzens_reaction/)")
    st.markdown("- [Arndt Eistert Reaction](https://www.organic-chemistry.org/namedreactions/arndt-eistert-synthesis.shtm/)")
    st.markdown("- [Hoffmann-Martius Rearrangement](https://en.wikipedia.org/wiki/Hofmann%E2%80%93Martius_rearrangement/)")
    st.markdown("### Named Reactions")
    st.markdown("- [Named Reactions in Organic Chemistry](https://organicchemistrydata.org/namedreactions/)")
    st.markdown("- [Comprehensive Named Reactions List](https://www.chem.ucla.edu/~harding/IGOC/N/named_reactions.html)")
    
    st.markdown("### Interactive Learning")
    st.markdown("- [Virtual Organic Chemistry Lab](https://www.labster.com/simulations/organic-chemistry/)")
    st.markdown("- [Chemical Structure Drawing](https://www.acdlabs.com/resources/freeware/#chemdrawdirect.perkinelmer.cloud/js/sample/index.html)")
    st.markdown("- [Chemical Structure Drawing](https://ccustomer.perkinelmer.com/#ncbi.nlm.nih.gov//edit3/index.html)")
    st.markdown("- [Chemical Structure Drawing](https//web.chemdoodle.com/)
    st.markdown("- [Chemical Structure Drawing](https//molview.org/)
    # Features explanation in sidebar
st.sidebar.title("Organic Synthesis Explorer")
st.sidebar.info("This app provides:")
st.sidebar.info("- Chemical equation balancing")
st.sidebar.info("- 50+ named organic reactions")
st.sidebar.info("- Synthesis planning resources")
st.sidebar.info("- Educational materials & links")

st.sidebar.divider()
st.sidebar.subheader("Featured Named Reactions")
st.sidebar.markdown("- Grignard Reaction")
st.sidebar.markdown("- Diels-Alder Reaction")
st.sidebar.markdown("- Wittig Reaction")
st.sidebar.markdown("- Friedel-Crafts Acylation")
st.sidebar.markdown("- Claisen Condensation")

st.sidebar.divider()
st.sidebar.subheader("Quick Links")
st.sidebar.markdown("[Organic Chemistry Portal](https://www.organic-chemistry.org/)")
st.sidebar.markdown("[Named Reactions Database](https://organicchemistrydata.org/namedreactions/)")
st.sidebar.markdown("[Virtual Chemistry Lab](https://www.labster.com/simulations/organic-chemistry/)")

st.sidebar.divider()
st.sidebar.caption("App developed by subramanianRamajayam &Powered by ChemPy for chemical balancing")
st.sidebar.caption("Educational resources from leading chemistry sources")
