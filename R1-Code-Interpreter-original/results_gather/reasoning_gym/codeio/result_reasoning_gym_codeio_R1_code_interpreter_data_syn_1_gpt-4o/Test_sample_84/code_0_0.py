def simulate_book_reading(book_title, author, nbr_pages, page_turns):
    # Start on the first page
    current_page = 1
    
    # Process each page turn
    for turn in page_turns:
        current_page += turn
        # Ensure the current page is within valid bounds
        if current_page < 1:
            current_page = 1
        elif current_page > nbr_pages:
            current_page = nbr_pages
    
    # Place a bookmark on the current page
    bookmarked_page = current_page
    
    # Return the result as a dictionary
    return {
        "current_page": current_page,
        "bookmarked_page": bookmarked_page
    }

# Given input
input_data = {'book_title': 'Book_70', 'author': 'Author_93', 'nbr_pages': 39, 'page_turns': [10]}

# Execute the function with the given input
result = simulate_book_reading(input_data['book_title'], input_data['author'], input_data['nbr_pages'], input_data['page_turns'])
print(result)