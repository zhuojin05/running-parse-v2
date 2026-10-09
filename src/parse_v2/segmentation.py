"""Domain segmentation and overlap resolution for ParSe v2."""

from typing import List, Literal, Tuple
from parse_v2.models import DomainRegion
from parse_v2.constants import DEFAULT_MIN_REGION_LEN, DEFAULT_PURITY_CUTOFF


def _scan_regions_for_label(
    labels: List[str],
    target_label: Literal["P", "D", "F"],
    min_len: int,
    purity_cutoff: float,
) -> List[Tuple[int, int]]:
    """Scan continuous regions with at least min_len and >= purity_cutoff of target_label.

    Returns list of 1-based (first_residue, last_residue) tuples.
    Matches Steve Whitten's Parse_v2.f logic exactly.
    """
    npep = len(labels)
    regions: List[Tuple[int, int]] = []
    
    i = 0  # 0-indexed in Python (corresponds to i=1 in Fortran)
    while (i + min_len - 1) < npep:
        count_target = 0
        count_w = 0
        region_found = False
        
        # Initial window of length min_len
        for j in range(i, i + min_len):
            count_w += 1
            if labels[j] == target_label:
                count_target += 1
                
        # Greedily extend to the right while purity >= purity_cutoff
        p_start = i
        p_end = i + min_len - 1
        
        while True:
            percent = count_target / count_w
            if percent >= purity_cutoff:
                region_found = True
                p_end = (i + count_w - 1)
                if p_end < (npep - 1):
                    # Peek next residue
                    count_w += 1
                    if labels[i + count_w - 1] == target_label:
                        count_target += 1
                    continue
                else:
                    break
            else:
                break
                
        if region_found:
            # 1-based indexing for residue numbers
            regions.append((p_start + 1, p_end + 1))
            i = p_end  # Next window starts at j in Fortran (0-indexed p_end)
        else:
            i += 1
            
    return regions


def segment_domains(
    labels: List[str],
    min_len: int = DEFAULT_MIN_REGION_LEN,
    purity_cutoff: float = DEFAULT_PURITY_CUTOFF,
) -> List[DomainRegion]:
    """Segment a sequence of residue labels into Folded, Disordered, and Phase-Separating domains.

    Includes overlap resolution matching Parse_v2.f.
    """
    p_regions = _scan_regions_for_label(labels, "P", min_len, purity_cutoff)
    d_regions = _scan_regions_for_label(labels, "D", min_len, purity_cutoff)
    f_regions = _scan_regions_for_label(labels, "F", min_len, purity_cutoff)

    # Combine all regions with their type: [type, start, end]
    all_domains = []
    for s, e in p_regions:
        all_domains.append(["P", s, e])
    for s, e in d_regions:
        all_domains.append(["D", s, e])
    for s, e in f_regions:
        all_domains.append(["F", s, e])

    if not all_domains:
        return []

    # Sort domains by start position
    all_domains.sort(key=lambda x: x[1])

    # Overlap resolution if more than 1 domain
    # Overlap split calculation: int((1.0 - purity_cutoff) * min_len / 2.0)
    split_offset = int((1.0 - purity_cutoff) * float(min_len) / 2.0)

    for i in range(len(all_domains) - 1):
        cur_type, cur_start, cur_end = all_domains[i]
        next_type, next_start, next_end = all_domains[i + 1]

        if cur_end >= next_start:
            # Domains overlap: truncate earlier domain and shift next domain
            new_cur_end = next_start + split_offset
            new_next_start = new_cur_end + 1
            all_domains[i][2] = new_cur_end
            all_domains[i + 1][1] = new_next_start

    # Build typed DomainRegion objects
    results = []
    for idx, (dtype, start, end) in enumerate(all_domains, 1):
        length = end - start + 1
        results.append(
            DomainRegion(
                region_id=idx,
                domain_type=dtype,
                first_residue=start,
                last_residue=end,
                length=length,
            )
        )

    return results
