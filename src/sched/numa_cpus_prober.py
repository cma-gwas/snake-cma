import os
import psutil
import numa


def get_numa_physical_cpus():
    cpu_affinities = os.sched_getaffinity(0)
    num_physical_cpus = psutil.cpu_count(logical=False)
    numa_node_cpus = numa.info.numa_hardware_info()['node_cpu_info']

    if len(cpu_affinities) >= num_physical_cpus * 2:
        pus_per_node = num_physical_cpus // numa.info.get_num_configured_nodes()
        node_pus = [numa_node_cpus[n][:pus_per_node] for n in numa_node_cpus.keys()]
        numa_pus = node_pus[0]
        for i in range(1, len(node_pus)):
            numa_pus.extend(node_pus[i])
        return numa_pus
    else:
        return cpu_affinities


if __name__ == "__main__":
    print(get_numa_physical_cpus())