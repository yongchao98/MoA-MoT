import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest single ACL entry (supernet) to match a list of networks.
    """
    # 1. The networks given in the problem.
    networks_str = ['172.20.96.0/19', '172.20.128.0/19']

    try:
        # 2. Convert the string representations into IPv4Network objects.
        networks = [ipaddress.ip_network(n) for n in networks_str]

        # 3. Calculate the smallest single network that contains all the given networks.
        # This is the supernet.
        summary_network = ipaddress.supernet(networks)

        # 4. Extract the network address and the wildcard mask.
        # The wildcard mask is the 'hostmask' in the ipaddress library.
        network_address = summary_network.network_address
        wildcard_mask = summary_network.hostmask

        # 5. Print the final result in the required format.
        # The numbers in the final equation are the network address and the wildcard mask.
        print(f"The smallest ACL entry that matches both networks is:")
        print(f"{network_address} {wildcard_mask}")

    except ValueError as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_summary_acl()