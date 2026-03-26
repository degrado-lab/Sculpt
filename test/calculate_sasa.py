import argparse
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

def calculate_residue_sasa(cif_file, chain_id, resid):
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', cif_file)
    except Exception as e:
        print(f"Error reading CIF file: {e}")
        return None
    
    # Calculate SASA for the whole structure
    # This computes SASA for all atoms/residues, since context matters
    sr = ShrakeRupley()
    sr.compute(structure, level="R") # Residue level computation
    
    # Extract SASA for the specific residue
    try:
        resid_int = int(resid)
        
        # Traverse to the requested chain and residue
        model = structure[0] # Assuming first model
        chain = model[chain_id]
        
        # Find the residue with the matching resseq (sequence identifier)
        target_residue = None
        for res in chain:
            # res.id is a tuple like (' ', 129, ' ') containing (hetero, resseq, icode)
            if res.id[1] == resid_int:
                target_residue = res
                break
                
        if target_residue is None:
            print(f"Error: Residue {resid} not found in chain {chain_id}.")
            return None
            
        # The ShrakeRupley compute() method adds a 'sasa' attribute to each residue
        sasa_value = target_residue.sasa
        print(f"SASA for residue {resid} in chain {chain_id} ({target_residue.get_resname()}): {sasa_value:.2f} Å²")
        return sasa_value
        
    except KeyError:
        print(f"Error: Chain {chain_id} not found in the structure.")
        return None
    except ValueError:
        print(f"Error: Invalid residue ID '{resid}'. Must be an integer.")
        return None
    except Exception as e:
        print(f"An error occurred while finding the residue: {e}")
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate SASA for a specific residue in a CIF file using Biopython.")
    parser.add_argument("cif_file", help="Path to the input CIF file")
    parser.add_argument("--chain", required=True, help="Chain ID (e.g., 'A')")
    parser.add_argument("--resid", required=True, help="Residue sequence ID (e.g., '129')")
    
    args = parser.parse_args()
    
    calculate_residue_sasa(args.cif_file, args.chain, args.resid)
