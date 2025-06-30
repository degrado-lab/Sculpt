### The following are the Scoring function classes, used to score a standardized structure file.
### Currently, only the geometry optimizer is implemented.
import tempfile
import ribbon

class SculptGeometricScoringFunction:
    """Geometric scoring function for evaluating molecular structures.
    
    This class wraps a user-defined function that takes a structure file as input
    and returns a geometric score based on structural features.
    """
    
    def __init__(self, function: callable):
        """Initialize the scoring function.
        
        Args:
            function: A callable that takes exactly one parameter (structure file)
                     and returns a numeric score.
                     
        Raises:
            ValueError: If the function doesn't take exactly one parameter.
        """
        self.function = function
        # verify function only has a single parameter:
        if function.__code__.co_argcount != 1:
            raise ValueError("Function must take exactly one parameter, an input structure file.")

    def score(self, structure_file):
        """Calculate the score for a structure file.
        
        Args:
            structure_file: Path to the structure file to score.
            
        Returns:
            The numeric score returned by the wrapped function.
        """
        # Call the function with the structure file
        return self.function(structure_file)

    def __call__(self, structure_file):
        """Allow the scoring function to be called directly.
        
        Args:
            structure_file: Path to the structure file to score.
            
        Returns:
            The numeric score returned by the wrapped function.
        """
        return self.score(structure_file)


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

    def score(self, structure_file):
        """Calculate the bond distance score for a structure.
        
        Args:
            structure_file: Path to the structure file to analyze.
            
        Returns:
            The difference between actual and target bond distance, optionally as absolute value.
            
        Raises:
            ValueError: If the distance file is empty or malformed.
        """
        # Calculate the score using ribbon
        outfile = tempfile.NamedTemporaryFile(delete=False, suffix='.dist')
        ribbon.CalculateDistance(
            pdb_file=structure_file,
            atom1=str(self.bond.atom_1),
            atom2=str(self.bond.atom_2),
            output_file=outfile.name,
        ).run()

        with open(outfile.name, 'r') as f:
            lines = f.readlines()
        if len(lines) < 1:
            raise ValueError("Distance file is empty or malformed.")
        distance = float(lines[0].strip())

        # To get the score, we get the difference between the target distance and the actual distance:
        score = distance - self.bond.target_distance

        if self.absolute:
            score = abs(score)

        return score

    def __call__(self, structure_file):
        """Allow the bond scorer to be called directly.
        
        Args:
            structure_file: Path to the structure file to analyze.
            
        Returns:
            The bond distance score.
        """
        return self.score(structure_file)

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

    def score(self, structure_file):
        """Calculate the bond angle score for a structure.
        
        Args:
            structure_file: Path to the structure file to analyze.
            
        Returns:
            The difference between actual and target bond angle, optionally as absolute value.
            
        Raises:
            ValueError: If the angle file is empty or malformed.
        """
        # Calculate the score using ribbon
        outfile = tempfile.NamedTemporaryFile(delete=False, suffix='.angle')
        ribbon.CalculateAngle(
            pdb_file=structure_file,
            atom1=str(self.angle.atom_1),
            atom2=str(self.angle.atom_2),
            atom3=str(self.angle.atom_3),
            output_file=outfile.name,
        ).run()

        with open(outfile.name, 'r') as f:
            lines = f.readlines()
        if len(lines) < 1:
            raise ValueError("Angle file is empty or malformed.")
        angle_value = float(lines[0].strip())

        # To get the score, we get the difference between the target angle and the actual angle:
        score = angle_value - self.angle.target_angle

        if self.absolute:
            score = abs(score)

        return score

    def __call__(self, structure_file):
        """Allow the angle scorer to be called directly.
        
        Args:
            structure_file: Path to the structure file to analyze.
            
        Returns:
            The bond angle score.
        """
        return self.score(structure_file)