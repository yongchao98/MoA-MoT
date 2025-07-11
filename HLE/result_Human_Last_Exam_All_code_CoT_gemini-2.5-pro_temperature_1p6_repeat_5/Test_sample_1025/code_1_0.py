def solve_reranking_task():
    """
    Analyzes two documents against two questions to determine their relative ranking.
    """
    
    # Storing the key information for clarity
    question1 = "What animal's coat color changed as a baby?"
    q1_doc1_evidence = "Young bulls changed their coat colour at a few months old..."
    q1_doc2_evidence = "young bulls changed their coat colour to black..."
    
    question2 = "What animal's coat color did not change as a baby?"
    q2_doc1_evidence = "Cows retained the reddish-brown colour."
    q2_doc2_evidence = "Cows retained the reddish-brown colour."

    print("Step 1: Analyze Question 1 relevance")
    print(f"  - Question: '{question1}'")
    print(f"  - Doc 1 Evidence: '{q1_doc1_evidence}'")
    print(f"  - Doc 2 Evidence: '{q1_doc2_evidence}'")
    print("  - Conclusion: Both documents state that the young bull's coat changed. They are equally relevant.")
    # Ranking for Q1 is equal.

    print("\nStep 2: Analyze Question 2 relevance")
    print(f"  - Question: '{question2}'")
    print(f"  - Doc 1 Evidence: '{q2_doc1_evidence}'")
    print(f"  - Doc 2 Evidence: '{q2_doc2_evidence}'")
    print("  - Conclusion: Both documents provide the identical sentence stating that the cow's coat color was retained. They are equally relevant.")
    # Ranking for Q2 is equal.
    
    print("\nStep 3: Final Conclusion")
    print("  - For both questions, the documents provide equally relevant information.")
    print("  - The correct option is 'Both are equally ranked'.")

    # Mapping the conclusion to the answer choices
    final_choice = "E"
    
    print(f"\nFinal Answer Choice: {final_choice}")

solve_reranking_task()