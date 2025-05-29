import re
from typing import List, Dict

def sort_classes(classes: List[str]) -> List[str]:
    # Regular expression to parse class definitions
    RE_CLASS = re.compile(r"class\s+(?P<class>\w+)(\((?P<parent>\w*)\))?")
    
    # Parse the class definitions to build a graph
    graph: Dict[str, List[str]] = {}
    indegree: Dict[str, int] = {}
    
    for c in classes:
        match = RE_CLASS.match(c)
        if match:
            class_name = match.group("class")
            parent_name = match.group("parent") or ""
            
            if class_name not in graph:
                graph[class_name] = []
            if parent_name and parent_name not in graph:
                graph[parent_name] = []
            
            if parent_name:
                graph[parent_name].append(class_name)
                indegree[class_name] = indegree.get(class_name, 0) + 1
            else:
                indegree[class_name] = indegree.get(class_name, 0)
    
    # Topological sort using Kahn's algorithm
    sorted_classes = []
    zero_indegree = [node for node in graph if indegree.get(node, 0) == 0]
    
    while zero_indegree:
        node = zero_indegree.pop()
        sorted_classes.append(node)
        
        for neighbor in graph[node]:
            indegree[neighbor] -= 1
            if indegree[neighbor] == 0:
                zero_indegree.append(neighbor)
    
    # Convert sorted class names back to class definitions
    class_definitions = {match.group("class"): c for c in classes for match in [RE_CLASS.match(c)] if match}
    return [class_definitions[class_name] for class_name in sorted_classes if class_name in class_definitions]

# Test the function with the given input
classes = ['class B', 'class A(B)', 'class E', 'class G']
sorted_classes = sort_classes(classes)
print(sorted_classes)