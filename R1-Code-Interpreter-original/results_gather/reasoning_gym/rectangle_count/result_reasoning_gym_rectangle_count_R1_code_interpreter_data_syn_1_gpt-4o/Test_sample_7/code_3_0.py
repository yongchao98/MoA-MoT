def count_rectangles(grid):
    # Convert the grid into a list of lines
    lines = grid.strip().split('\n')
    
    # Initialize variables
    rectangles = 0
    in_rectangle = False
    overlap_detected = False
    
    # Iterate over each line in the grid
    for line in lines:
        # Check if the line contains a rectangle boundary
        if '#' in line or '█' in line:
            if not in_rectangle:
                # We found the start of a new rectangle
                rectangles += 1
                in_rectangle = True
            if '█' in line and not overlap_detected:
                # We found an overlap, count the overlapping rectangle
                rectangles += 1
                overlap_detected = True
        else:
            # No rectangle boundary in this line, reset the flags
            in_rectangle = False
            overlap_detected = False
    
    return rectangles

# Define the grid
grid = """
                                                                                
                                                                                
                                                                                
                                 ###########################                    
                                 #                         #                    
                                 #                         #                    
                                 #                         #                    
                                 #                         #                    
                                 #                         #                    
   ##############################█#########                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #        #                #                    
   #                             #########█#################                    
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                                      #                                     
   #                            ##########█###################################  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #                            #         #                                  #  
   #############################█##########                                  #  
                                ##############################################  
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                            ############################                        
                            #                          #                        
                            ############################                        
                                                                                
                                                                                
                                                                                
                                                                                
"""

# Count the rectangles
rectangles_count = count_rectangles(grid)
print(rectangles_count)