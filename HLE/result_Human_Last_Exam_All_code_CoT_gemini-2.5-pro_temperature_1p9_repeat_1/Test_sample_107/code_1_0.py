import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest ACL entry to summarize two given IP networks.
    """
    try:
        # Define the two networks to be summarized
        net1_str = '172.20.96.0/19'
        net2_str = '172.20.128.0/19'

        # Create ip_network objects for robust calculations
        network1 = ipaddress.ip_network(net1_str)
        network2 = ipaddress.ip_network(net2_str)

        # collapse_addresses finds the minimal set of networks covering the inputs
        # In this case, it will return a single summary network.
        summary_list = list(ipaddress.collapse_addresses([network1, network2]))
        
        if not summary_list:
            print("Could not summarize the given networks.")
            return

        # The result is the first network in the summarized list
        summary_net = summary_list[0]

        # The network address for the ACL entry
        network_address = summary_net.network_address

        # The wildcard mask is the inverse of the netmask, 
        # which is directly available as the hostmask in the ipaddress library.
        wildcard_mask = summary_net.hostmask

        # Print the final ACL entry with each number of the address and mask
        print("The smallest appropriate IP access control list entry is:")
        print(f"{network_address} {wildcard_mask}")

    except ValueError as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_summary_acl()