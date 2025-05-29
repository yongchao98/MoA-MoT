import json

def feedingTime(classes):
    N = len(classes)
    colorArr = [0] * N
    
    # Creating the graph out of the classes
    graph = [[0 for _ in range(N)] for _ in range(N)]  # 2D array for graph
    classes = [set(json.loads(c)) for c in classes]
    
    # Iterate through set and check where there are intersections
    for i in range(N-1):
        for j in range(i+1, N):
            if classes[i].intersection(classes[j]):
                graph[i][j] = 1
                graph[j][i] = 1
    
    def isColorable(curr, colorArr, color):
        for nei in range(N):
            if graph[curr][nei] == 1 and colorArr[nei] == color:
                return False  # A neighbor shares the color, not colorable with this color
        return True
    
    def coloringGraph(curr, k):
        if curr == N:
            return True  # It has recursed till past the last node and has passed the colorability check
        
        for color in range(1, k+1):  # See if k coloring works by trial and error of different colors
            if isColorable(curr, colorArr, color):
                colorArr[curr] = color  # Colored the current node
                
                # Now check neighbors recursively
                if coloringGraph(curr+1, k):
                    return True  # Get here if we iterate through all neighbors and curr = N at the top, then it passed coloring test
                
                # If above test fails, then we know that the neighbors are not k colorable, reset to try again in future
                colorArr[curr] = 0
    
        return False  # We haven't returned True until this point so the coloring range does not work with k, return False
                
        
    for k in range(1, 6):
        if coloringGraph(0, k):
            return k  # Found a k coloring that works
    
    return -1  # No coloring found within 5 colorable, so return -1

# Main function
def main_solution(classes):
    # Call the original function
    result = feedingTime(classes)
    
    # Convert the result to JSON serializable output
    return result

# Given input
input_data = {'classes': ['["zebra", "penguin", "monkey", "giraffe", "tiger", "elephant", "snake", "bear"]', 
                          '["bear", "kangaroo", "penguin", "monkey"]', 
                          '["giraffe", "penguin", "zebra", "monkey", "elephant", "bear", "kangaroo"]', 
                          '["elephant"]', 
                          '["snake", "kangaroo", "penguin", "lion"]', 
                          '["monkey", "tiger", "giraffe", "lion", "bear", "zebra", "snake"]']}

# Execute the solution
print(main_solution(input_data['classes']))