from pathlib import Path
from ppg_bp_database.datasets import bp_graphene_utils


if __name__ == "__main__":
    bp_graphene_utils.convert(Path("data/raw/bpgraphene"), Path("data/structured/bpgraphene"))
