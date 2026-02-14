#from sculpt.id.word_list import LEFT, RIGHT
from sculpt.id.word_list import LEFT, RIGHT
import random
import uuid

def generate_id(suffix=True):
    left = random.choice(LEFT)
    right = random.choice(RIGHT)
    name = f"{left}_{right}"
    if suffix:
        name += f"_{str(uuid.uuid4())[:8]}"
    return name