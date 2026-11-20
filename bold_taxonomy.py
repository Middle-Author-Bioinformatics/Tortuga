"""Normalize missing taxonomy values in current and previously prepared BOLD records."""

MISSING = {'', 'none', 'null', 'nan', 'na', 'n/a', '-', 'not provided'}


def taxon_name(value):
    text = str(value).strip() if value is not None else ''
    return '' if text.casefold() in MISSING else text


def normalize_taxonomy(taxon, lineage):
    """Recover the deepest observed name without inventing an expected identity."""
    names = [name for part in str(lineage or '').split(';') if (name := taxon_name(part))]
    return taxon_name(taxon) or (names[-1] if names else ''), ';'.join(names)
