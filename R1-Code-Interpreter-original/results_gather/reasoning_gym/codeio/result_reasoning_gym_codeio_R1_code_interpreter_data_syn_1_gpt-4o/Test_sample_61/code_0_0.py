import re
from collections import defaultdict

def merge_intervals(sorted_by_lower_bound):
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher[0] <= (lower[1] + 1):
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)
            else:
                merged.append(higher)
    return merged

def mapped_intervals(pos, cigar, md):
    if pos == "0":
        return []
    list_mapped_intervals = []
    current_pos = int(pos)
    
    if md != "":
        p = re.compile('\d+\D*')
        list_md = p.findall(md)
        for x in range(len(list_md)):
            p_num = re.compile('\d+')
            dist_match = int(p_num.findall(list_md[x])[0])
            if "^" in list_md[x]:
                dist_skip = len(list_md[x]) - len(str(dist_match)) - 1
            else:
                dist_skip = len(list_md[x]) - len(str(dist_match))
            if dist_match != 0:
                list_mapped_intervals.append((current_pos, current_pos + dist_match - 1))
            current_pos = current_pos + dist_match + dist_skip
    else:
        p = re.compile('\d+\D')
        list_cigar_ope = p.findall(cigar)
        for x in range(len(list_cigar_ope)):
            distance = int(list_cigar_ope[x][:-1])
            operation = list_cigar_ope[x][-1]
            if operation == "M":
                list_mapped_intervals.append((current_pos, current_pos + distance - 1))
                current_pos = current_pos + distance
            elif operation == "D":
                current_pos = current_pos + distance
    return list_mapped_intervals

def main_solution(pos, cigar, md):
    intervals = mapped_intervals(pos, cigar, md)
    sorted_intervals = sorted(intervals, key=lambda tup: tup[0])
    merged_intervals = merge_intervals(sorted_intervals)
    return [list(interval) for interval in merged_intervals]

# Constructed inputs
pos = "943"
cigar = "10M1D6M1D10M1D10M"
md = "10^A6^C10^G10"

# Verify the output
print(main_solution(pos, cigar, md))