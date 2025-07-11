import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest appropriate IP access control list entry (network address and wildcard mask)
    to match hosts on a given list of networks.
    """
    # List of networks in CIDR format
    network_strings = ['172.20.96.0/19', '172.20.128.0/19']

    # Convert the string representations into ip_network objects
    networks = [ipaddress.ip_network(net_str) for net_str in network_strings]

    # The collapse_addresses function finds the minimal set of networks
    # that cover all the addresses in the input list.
    # Since the input networks are adjacent or can be summarized,
    # this will give us the smallest possible supernet(s).
    summary_networks = list(ipaddress.collapse_addresses(networks))

    # For this specific case, we expect a single summary network.
    if len(summary_networks) == 1:
        supernet = summary_networks[0]
        
        # The network address of the summary network
        network_address = supernet.network_address
        
        # The hostmask is the equivalent of a wildcard mask in ACLs
        wildcard_mask = supernet.hostmask
        
        # Print the final result: network address followed by the wildcard mask
        print(f"{network_address} {wildcard_mask}")
    else:
        # Handle cases where summarization might result in multiple networks
        print("Could not find a single summary network. Resulting networks:")
        for net in summary_networks:
            print(f"{net.network_address} {net.hostmask}")

# Execute the function
calculate_summary_acl()