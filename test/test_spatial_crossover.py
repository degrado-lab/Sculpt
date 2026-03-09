
import sys
import os
from pathlib import Path

# Add Sculpt to sys.path
sys.path.append(str(Path('/home/freitas/workspace/sculpt_public_release/Sculpt')))

from sculpt.design import Design
from sculpt.ga.crossover import SpatialCrossover

def test_spatial_crossover():
    # Paths to example structure
    cif_path = Path('/home/freitas/workspace/sculpt_public_release/Sculpt/run_data/data/5RGA_fixed/5RGA_KEMP1_fixed.cif')
    
    if not cif_path.exists():
        print(f"Skipping test: {cif_path} not found.")
        return

    # Create two designs
    d1 = Design(name="parent1")
    d1.load_structure(cif_path)
    d1.load_sequence_from_structure()
    
    d2 = Design(name="parent2")
    d2.load_structure(cif_path)
    d2.load_sequence_from_structure()
    
    # Initialize SpatialCrossover
    sc = SpatialCrossover(crossover_rate=1.0)
    
    # Run crossover
    population = [d1, d2]
    offspring = sc.crossover(population, unique_naming=False)
    
    print(f"Offspring count: {len(offspring)}")
    for d in offspring:
        print(f"Design: {d.name}")
        if d.sequence:
            print(f"Sequence length: {len(d.sequence.split('\\n')[-1])}")
            print(f"History: {d.history}")
        else:
            print("No sequence generated.")
            
    assert len(offspring) == 2
    assert "Spatial crossover" in offspring[0].history[-1]

if __name__ == "__main__":
    try:
        test_spatial_crossover()
        print("Test passed!")
    except Exception as e:
        print(f"Test failed: {e}")
        import traceback
        traceback.print_exc()
