### The following are the Scoring function classes, used to score a standardized structure file.
### Currently, only the geometry optimizer is implemented.
import tempfile
import ribbon

from sculpt.geometry import CustomBond, CustomAngle, CustomTorsion
from sculpt.design import Design

class SculptGeometricScoringFunction:
    """Geometric scoring function for evaluating molecular structures.
    
    This class wraps a user-defined function that takes a structure file as input
    and returns a geometric score based on structural features.
    """
    
    def __init__(self, function: callable, scorers: dict = None, aggregator: callable = None):
        """Initialize the scoring function.
        
        Args:
            function: A callable that takes exactly one parameter (Design object)
                     and returns a numeric score.
            scorers: A dictionary of scorers to be used for scoring. (in the case of multi-state design)
            aggregator: A callable that takes a dictionary of scores and returns a single score. (in the case of multi-state design)
                     
        Raises:
            ValueError: If the function doesn't take exactly one parameter.
        """
        self.function = function
        self.scorers = scorers
        self.aggregator = aggregator

        # verify function only has a single parameter:
        if function.__code__.co_argcount != 1:
            raise ValueError("Function must take exactly one parameter, an input Design.")

    def score(self, design: Design):
        """Calculate the score for a structure file.
        
        Args:
            structure_file: Path to the structure file to score.
            
        Returns:
            The numeric score returned by the wrapped function.
        """
        if self.scorers and self.aggregator:
            scorer_results = {}
            for name, scorer in self.scorers.items():
                scorer_results[name] = scorer(design)
            return self.aggregator(scorer_results)

        elif self.function:
            # Call the function with the structure file
            return self.function(design)
        else:
            raise ValueError("No valid scoring mechanisms (function or scorers+aggregator) were provided.")

    def __call__(self, design: Design):
        """Allow the scoring function to be called directly.
        
        Args:
            structure_file: Path to the structure file to score.
            
        Returns:
            The numeric score returned by the wrapped function.
        """
        return self.score(design)


### The following classes define how we calculate a score for each design.
### For example, ScoreBond calculates the score for each bond in a design, using a provided CustomBond object.

class ScoreBond:
    """Score calculator for bond distances in molecular structures.
    
    This class calculates scores based on the deviation from target bond distances
    using the CustomBond restraints.
    """
    
    def __init__(self, bond: CustomBond, absolute=True):
        """Initialize the bond scorer.
        
        Args:
            bond: CustomBond object containing the bond definition and target distance.
            absolute: If True, return the absolute value of the score. Defaults to True.
        """
        self.bond = bond
        self.absolute = absolute # Get the absolute value of the score if True

    def score(self, design):
        """Calculate the bond distance score for a structure.
        
        Args:
            design: The structure object to analyze.

        Returns:
            The difference between actual and target bond distance, optionally as absolute value.
            
        Raises:
            ValueError: If the distance file is empty or malformed.
        """

        with tempfile.NamedTemporaryFile(delete=False, suffix='.dist') as outfile, \
             design.temp_structure_file() as structure_file:

            # Calculate the score using ribbon
            ribbon.CalculateDistance(
                pdb_file=structure_file,
                atom1=str(self.bond.atom_1),
                atom2=str(self.bond.atom_2),
                output_file=outfile.name,
            ).run()


            try:
                with open(outfile.name, 'r') as f:
                    lines = f.readlines()
                    print(lines)
                if len(lines) < 1:
                    raise ValueError("Distance file is empty or malformed.")
                distance = float(lines[0].strip())
            except Exception as e:
                distance = float('inf')
            

            # To get the score, we get the difference between the target distance and the actual distance:
            score = distance - self.bond.target_distance

            if self.absolute:
                score = abs(score)

            ### DEBUG
            print(f"Distance value: {distance}")
            print(f"Distance target: {self.bond.target_distance}")
            print(f"Distance score: {score}")

        return score

    def __call__(self, design):
        """Allow the bond scorer to be called directly.
        
        Args:
            structure_file: Path to the structure file to analyze.
            
        Returns:
            The bond distance score.
        """
        return self.score(design)

### N.B. Does this work with Ribbon yet?
class ScoreAngle:
    """Score calculator for bond angles in molecular structures.
    
    This class calculates scores based on the deviation from target bond angles
    using the CustomAngle restraints.
    """
    
    def __init__(self, angle: CustomAngle, absolute=True):
        """Initialize the angle scorer.
        
        Args:
            angle: CustomAngle object containing the angle definition and target angle.
            absolute: If True, return the absolute value of the score. Defaults to True.
        """
        self.angle = angle
        self.absolute = absolute  # Get the absolute value of the score if True

    def score(self, design: Design):
        """Calculate the bond angle score for a structure.
        
        Args:
            design: The Design object to analyze.

        Returns:
            The bond angle score.
        """
        with tempfile.NamedTemporaryFile(delete=False, suffix='.angle') as outfile, \
                design.temp_structure_file() as structure_file:
            # Calculate the score using ribbon
            ribbon.CalculateAngle(
                pdb_file=structure_file,
                atom1=str(self.angle.atom_1),
                atom2=str(self.angle.atom_2),
                atom3=str(self.angle.atom_3),
                output_file=outfile.name,
            ).run()

            try:
                with open(outfile.name, 'r') as f:
                    lines = f.readlines()
                if len(lines) < 1:
                    raise ValueError("Angle file is empty or malformed.")
                angle_value = float(lines[0].strip())
            except Exception as e:
                angle_value = float('inf')

            # To get the score, we get the difference between the target angle and the actual angle:
            score = angle_value - self.angle.target_angle

            if self.absolute:
                score = abs(score)

            ### DEBUG
            print(f"Angle value: {angle_value}")
            print(f"Angle target: {self.angle.target_angle}")
            print(f"Angle score: {score}")

        return score

    def __call__(self, design: Design):
        """Allow the angle scorer to be called directly.
        
        Args:
            design: The Design object to analyze.
        Returns:
            The bond angle score.
        """
        return self.score(design)

class ScoreDihedral:
    """Score calculator for dihedral angles in molecular structures.
    
    This class calculates scores based on the deviation from target dihedral angles
    using the CustomTorsion restraints.
    """
    
    def __init__(self, torsion: CustomTorsion, absolute=True):
        """Initialize the dihedral scorer.
        
        Args:
            torsion: CustomTorsion object containing the torsion definition and target angle.
            absolute: If True, return the absolute value of the score. Defaults to True.
        """
        self.torsion = torsion
        self.absolute = absolute  # Get the absolute value of the score if True
        self.periodicity = torsion.periodicity 

    def score(self, design: Design):
        """Calculate the dihedral angle score for a structure.
        
        Args:
            design: The Design object to analyze.

        Returns:
            The dihedral angle score.
        """
        with tempfile.NamedTemporaryFile(delete=False, suffix='.dihedral') as outfile, \
                design.temp_structure_file() as structure_file:
            # Calculate the score using ribbon
            ribbon.CalculateDihedral(
                pdb_file=structure_file,
                atom1=str(self.torsion.atom_1),
                atom2=str(self.torsion.atom_2),
                atom3=str(self.torsion.atom_3),
                atom4=str(self.torsion.atom_4),
                output_file=outfile.name,
            ).run()

            try:
                with open(outfile.name, 'r') as f:
                    lines = f.readlines()
                if len(lines) < 1:
                    raise ValueError("Dihedral file is empty or malformed.")
                dihedral_value = float(lines[0].strip())
            except Exception as e:
                dihedral_value = float('inf')

            # To get the score, we get the difference between the target angle and the actual angle:
            score = dihedral_value - self.torsion.target_angle

            # Apply the periodicity correction if needed
            # Not yet implemented

            if self.absolute:
                score = abs(score)

            ### DEBUG
            print(f"Dihedral value: {dihedral_value}")
            print(f"Dihedral target: {self.torsion.target_angle}")
            print(f"Dihedral score: {score}")

        return score

    def __call__(self, design: Design):
        """Allow the dihedral scorer to be called directly.

        Args:
            design: The Design object to analyze.
        Returns:
            The dihedral angle score.
        """
        return self.score(design)