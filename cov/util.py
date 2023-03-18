from pathlib import Path
from collections import namedtuple


def get_data_path():
    """
    Get path to data folder in base directory
    """
    return Path(__file__).parents[1] / "data"


def get_data_folders(data_path=None):
    """
    Generate named tuple for the path to processed, raw, and tables folder.
    """
    if data_path is None:
        data_path = get_data_path()
    fields = ["processed", "raw", "tables"]
    dir = namedtuple("dir", fields)
    data_fold = dir(*[data_path / f for f in fields])
    return data_fold


def get_raw_folders(dataset=None):
    """
    Args:
        dataset (str, optional): One of ['human_OE', 'mouse_OE', 'mouse_OB', 'mouse_HBC'] for the respective datasets.
        If none, returns all paths.

    Returns:
        Pathlib.Path: Path to raw data folder for each dataset
    """
    data_fold = get_data_folders()
    raw_fold = data_fold.raw
    mapping = {
        "human_OE": "GSE139522_RAW",
        "mouse_OE": "GSE151346_MOE",
        "mouse_OB": "GSE148360_OB",
        "mouse_HBC": "GSE157068_K5",
    }
    raw_folds = {k: raw_fold / v for k, v in mapping.items()}
    if dataset is None:
        return raw_folds
    elif dataset in raw_folds:
        return raw_folds[dataset]
    else:
        raise ValueError(f"Dataset {dataset} not one of {mapping.keys()}")
