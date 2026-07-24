"""Versioned constants and controlled vocabularies."""

PIPELINE_NAME = "biosensia-candidate-pocket-library"
PIPELINE_VERSION = "0.2.0"
MANIFEST_SCHEMA_VERSION = "2.0.0"
SIDECAR_SCHEMA_VERSION = "2.0.0"
EXTRACTION_SCHEMA_VERSION = "1"
LIGAND_PARSER_SCHEMA_VERSION = "1"
HASH_SCHEMA_VERSION = "1"
PICKLE_PROTOCOL = 4

CANONICAL_AMINO_ACIDS = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL",
    }
)
DEFAULT_MODIFIED_AMINO_ACIDS = ("MSE",)
WATER_RESIDUES = frozenset({"HOH", "WAT", "DOD"})
YEAR_RANGES = ("1981-2000", "2001-2010", "2011-2019")
EXPECTED_COMPLEX_SUFFIXES = (
    "_ligand.sdf", "_ligand.mol2", "_pocket.pdb", "_protein.pdb"
)
INDEX_FILES = (
    "INDEX_general_NL.2020R1.lst",
    "INDEX_general_PL.2020R1.lst",
    "INDEX_general_PN.2020R1.lst",
    "INDEX_general_PP.2020R1.lst",
    "README",
)
SUPPORTED_GEOMETRY_TIERS = frozenset({"A", "B", "C"})
ISSUE_SEVERITIES = frozenset({"info", "warning", "error", "fatal"})
