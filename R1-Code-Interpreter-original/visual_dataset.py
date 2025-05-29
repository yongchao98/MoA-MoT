import json
import os


def display_fixed_examples(json_path, indices=[0, 1, 2]):
    """
    Load a JSON dataset, and display examples for the given list of indices.
    """
    if not os.path.isfile(json_path):
        print(f"File not found: {json_path}")
        return

    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    print(f"Loaded {len(data)} items from: {json_path}\n")

    if not data:
        print("The dataset is empty.")
        return

    for i in indices:
        if i >= len(data):
            print(f"Index {i} is out of range for dataset length {len(data)}.\n")
            continue

        example = data[i]
        print(f"===== EXAMPLE {i} =====")

        # Display the conversation history
        history = example.get("history", [])
        if history:
            print("Conversation History:")
            for turn_num, turn in enumerate(history, start=1):
                user_msg, assistant_msg = turn
                print(f"  [Turn {turn_num}] USER:\n    {user_msg}")
                print(f"  [Turn {turn_num}] ASSISTANT:\n    {assistant_msg}")
        else:
            print("No conversation history found for this item.")

        # Display the final instruction
        print("\nFinal Instruction:")
        print(f"  {example.get('instruction', '')}")

        # Display the final output
        print("\nFinal Output:")
        print(f"  {example.get('output', '')}")

        print("\n" + "=" * 25 + "\n")


if __name__ == '__main__':
    # Update this path to point to your JSON file
    json_file_path = '/home/ycchen/R1-Code-Interpreter/results_gather/R1_CI_gather_dataset_round2.json'

    # By default, shows items at indices 0, 1, 2
    display_fixed_examples(json_file_path, indices=[0,1,2])
