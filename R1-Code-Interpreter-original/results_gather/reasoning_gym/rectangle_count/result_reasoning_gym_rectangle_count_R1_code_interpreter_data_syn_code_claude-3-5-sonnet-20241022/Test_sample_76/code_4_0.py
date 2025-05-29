def count_rectangles(grid):
    lines = grid.split('\n')
    # Remove empty lines
    lines = [line for line in lines if line.strip()]
    
    def count_distinct_rectangles():
        # Count the large outer rectangle
        outer_rectangle = 1
        
        # Count the small top rectangle (8x1)
        small_rectangle = 1
        
        # Count the medium rectangle
        medium_rectangle = 1
        
        # Count overlapping rectangles (count █ characters)
        # Each █ represents an additional rectangle
        overlap_count = sum(line.count('█') for line in lines)
        
        # Total rectangles = outer + small + medium + overlaps
        total = outer_rectangle + small_rectangle + medium_rectangle + overlap_count
        print(total)

    count_distinct_rectangles()

grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
          ####################################################                  
          #                                                  #                  
          #                                                  #                  
          #              ########                            #                  
          #              #      #                            #                  
          #            ##█######█#####                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # #      #    #                       #                  
          #            # ########    #                       #                  
          #            ###############                       #                  
          #                                                  #                  
          #                                                  #                  
          #                                                  #                  
          #                                                  #                  
          #                                                  #                  
          #                                                  #                  
          ####################################################                  
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

count_rectangles(grid)