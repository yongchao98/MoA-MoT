# Initialize the grid as a string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
            ########################                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            #                      #                                            
            ########################                                            
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

# Convert the grid into a list of lines
lines = grid.split('\n')

# Count rectangles
count = 0

# Look for complete rectangles
for i in range(len(lines)):
    for j in range(len(lines[i])):
        if lines[i][j] == '#':  # Found potential top-left corner
            # Look for complete rectangle
            width = 0
            height = 0
            
            # Find width
            for k in range(j, len(lines[i])):
                if lines[i][k] == '#':
                    width += 1
                else:
                    break
            
            # Find height
            for k in range(i, len(lines)):
                if k < len(lines) and j < len(lines[k]) and lines[k][j] == '#':
                    height += 1
                else:
                    break
            
            # Verify if it's a complete rectangle
            is_rectangle = True
            for y in range(i, i + height):
                for x in range(j, j + width):
                    if y >= len(lines) or x >= len(lines[y]) or lines[y][x] != '#':
                        is_rectangle = False
                        break
                if not is_rectangle:
                    break
            
            # If it's a complete rectangle and we're at its top-left corner
            if (is_rectangle and 
                (i == 0 or j == 0 or lines[i-1][j] != '#') and 
                (i == 0 or j == 0 or lines[i][j-1] != '#')):
                count += 1

print(count)