class Link:
    empty = ()
    
    def __init__(self, first, rest=empty):
        assert rest is Link.empty or isinstance(rest, Link)
        self.first = first
        self.rest = rest
        
    def to_list(self):
        """Convert the linked list to a Python list."""
        if self.rest is Link.empty:
            return [self.first]
        else:
            return [self.first] + self.rest.to_list()

def sorted_list_to_linked_list(sorted_list):
    """Helper function to convert a sorted list to a linked list."""
    if not sorted_list:
        return Link.empty
    else:
        return Link(sorted_list[0], sorted_list_to_linked_list(sorted_list[1:]))

def main_solution(elements):
    # Sort the input list
    sorted_elements = sorted(elements)
    
    # Convert the sorted list to a linked list
    sorted_linked_list = sorted_list_to_linked_list(sorted_elements)
    
    # Convert the linked list back to a Python list and return
    return sorted_linked_list.to_list()

# Given input
elements = [79, 17, 15, -66, -18, -73, -61, 72, 72, 53]

# Execute the main solution function and print the result
print(main_solution(elements))