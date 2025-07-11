def find_author():
    """
    This function provides the answer to the user's question about a classical quote.
    """
    author = "Quintilian"
    work = "Institutio Oratoria (Institutes of Oratory), Book 2, Chapter 2"
    quote = "prope soli iam in scholis erunt relicti"
    explanation = (
        "The Roman rhetorician Quintilian used a version of this quote in his work, "
        f"'{work}'.\n\n"
        "He was discussing the challenges for teachers of oratory. He argued that "
        "if professors do not offer profound and valuable insights beyond a basic curriculum, "
        "they risk becoming irrelevant to their students. In his words, they would be "
        "'left almost alone in their schools' ('prope soli iam in scholis erunt relicti'). "
        "This perfectly matches the context of your question."
    )

    print(f"The quote is from the classical author: {author}")
    print("\n--- Details ---")
    print(explanation)

if __name__ == '__main__':
    find_author()