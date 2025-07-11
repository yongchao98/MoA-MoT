import collections

def get_genji_ko_partition(chapter_number):
    """
    Returns the Genji-kō partition for a given chapter.
    This function stores the data for chapter 39.
    """
    # Genji-kō patterns data (only for chapter 39 as required)
    # The pattern for Chapter 39, "Yūgiri", connects {1, 2, 3} and {4, 5}.
    patterns = {
        39: [[1, 2, 3], [4, 5]]
    }
    return patterns.get(chapter_number)

def format_partition(partition):
    """
    Formats the partition into a string representation of a sorted set of sets.
    e.g., [[1,3,4], [2], [5]] -> "{{1,3,4},{2},{5}}"
    """
    # Sort elements within each subset
    sorted_subsets = [sorted(subset) for subset in partition]
    
    # Sort the subsets themselves based on their first element
    # Python's default list sort does this automatically.
    sorted_subsets.sort()
    
    # Format each subset into a string like "{1,2,3}"
    subset_strings = []
    for subset in sorted_subsets:
        # Create the string for each number in the subset
        number_strings = [str(num) for num in subset]
        subset_strings.append(f"{{{','.join(number_strings)}}}")
        
    # Join the formatted subset strings and wrap with outer braces
    return f"{{{','.join(subset_strings)}}}"

# Chapter in question
chapter = 39

# Get the partition for the chapter
partition = get_genji_ko_partition(chapter)

# Format and print the result
if partition:
    formatted_result = format_partition(partition)
    print(f"The incense pattern partition for chapter {chapter} corresponds to: {formatted_result}")
else:
    print(f"Partition data for chapter {chapter} is not available.")
