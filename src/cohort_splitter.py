import time
import numpy as np
import multiprocessing as mp
import argparse
import psutil
import pandas as pd


def _write_chunk(arr, start_idx, end_idx, output_path):
    #np.savetxt(output_path, np.sort(arr[start_idx:end_idx], axis=0), delimiter=' ', fmt='%d')
    np.savetxt(output_path, np.sort(arr[start_idx:end_idx], axis=0), delimiter=' ', fmt='%s')
    return output_path


def _write_chunk_bt(arr, start_idx, end_idx, case_arr, case_start_idx, case_end_idx, output_path):
    sub_arr = np.concatenate([arr[start_idx:end_idx], case_arr[case_start_idx:case_end_idx]])

    #np.savetxt(output_path, np.sort(sub_arr, axis=0), delimiter=' ', fmt='%d')
    np.savetxt(output_path, np.sort(sub_arr, axis=0), delimiter=' ', fmt='%s')
    return output_path


def write_splits(samples_arr, n_splits, output_prefix, n_threads=1):
    split_size, rem = divmod(len(samples_arr), n_splits)

    start_time = time.perf_counter()
    with mp.Pool(n_threads) as pool:
        async_refs = []
        for s in range(n_splits):
            start, end = s * split_size, (s+1) * split_size
            if rem > 0 and s == n_splits - 1:
                end += rem

            async_refs.append(pool.apply_async(_write_chunk, args=(samples_arr, start, end,
                                                                   output_prefix + ".split.{}".format(s+1))))

        [print(ref.get()) for ref in async_refs]

    print(" Splitting finished ({:.2f} s)".format(time.perf_counter() - start_time))


def write_splits_bt(control_arr, case_arr, n_splits, output_prefix, overlap_cases=False, n_threads=1):
    split_size, rem = divmod(len(control_arr), n_splits)
    if not overlap_cases:
        case_split_size, case_rem = divmod(len(case_arr), n_splits)

    start_time = time.perf_counter()
    with mp.Pool(n_threads) as pool:
        async_refs = []
        for s in range(n_splits):
            start, end = s * split_size, (s+1) * split_size
            if rem > 0 and s == n_splits - 1:
                end += rem

            if not overlap_cases:
                case_start, case_end = s * case_split_size, (s+1) * case_split_size
                if case_rem > 0 and s == n_splits - 1:
                    case_end += case_rem
            else:
                case_start, case_end = 0, len(case_arr)

            async_refs.append(pool.apply_async(_write_chunk_bt, args=(control_arr, start, end, case_arr, case_start,
                                                                      case_end, output_prefix + ".split.{}".format(s+1))))

        [print(ref.get()) for ref in async_refs]

    print(" Splitting finished ({:.2f} s)".format(time.perf_counter() - start_time))


def write_bisect(samples_arr, output_prefix, bisect_ratio, n_threads=1):
    cut_offset = int(bisect_ratio * len(samples_arr))

    start_time = time.perf_counter()
    with mp.Pool(n_threads) as pool:
        async_refs = list()

        async_refs.append(pool.apply_async(_write_chunk, args=(samples_arr, 0, cut_offset, output_prefix + ".split.1")))
        async_refs.append(pool.apply_async(_write_chunk, args=(samples_arr, cut_offset, len(samples_arr),
                                                               output_prefix + ".split.2")))

        [print(ref.get()) for ref in async_refs]

    print(" Splitting finished ({:.2f} s)".format(time.perf_counter() - start_time))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample_file", type=str, help="Sample file to be split", required=True)
    parser.add_argument("-o", "--output_prefix", type=str, help="Prefix for output path", required=True)
    parser.add_argument('-r', '--randomize', action='store_true')  # on/off flag
    parser.add_argument("-n", "--num_splits", type=int, help="Number of splits", default=2)
    parser.add_argument("--rand_seed", type=int, help="Random seed", default=None)
    parser.add_argument("--threads", type=int, help="Number of parallel threads", default=None)
    parser.add_argument("--bisect-ratio", type=float, help="Split ratio for bisect splitting")
    parser.add_argument("--sorted-id-list", type=str, help="Sorted list of individuals")
    parser.add_argument("--binary-traits", type=str, help="Phenotype file used to balance case-control")
    parser.add_argument("--overlap-case-samples", action='store_true', help="Add case samples to all splits")
    args = parser.parse_args()

    sample_file_path = args.sample_file
    output_prefix = args.output_prefix
    num_splits = args.num_splits

    if args.rand_seed:
        rng = np.random.default_rng(args.rand_seed)
    else:
        rng = np.random.default_rng()

    # if args.randomize:
    #     if args.rand_seed:
    #         rng = np.random.default_rng(args.rand_seed)
    #     else:
    #         rng = np.random.default_rng()
    # else:
    #     print(" --rand_seed is ignored because --randomize is False")

    if args.bisect_ratio:
        assert 0 < args.bisect_ratio < 1.0, " --bisect_ratio must be in the range [0, 1]"

        print(" --num_splits is ignored because --bisect-ratio is set")
        num_splits = 2

    if args.threads:
        num_threads = args.threads
    else:
        num_threads = psutil.cpu_count(logical=False)

    #sample_id_arr = np.genfromtxt(sample_file_path, usecols=(0, 1))
    #print(sample_id_arr)

    if args.sorted_id_list:
        #sorted_id_arr = np.genfromtxt(args.sorted_id_list, usecols=(0, 1))
        sorted_id_df = pd.read_csv(args.sorted_id_list, delim_whitespace=True, header=None, usecols=[0],
                                   dtype='string')
        sorted_id_df["IID"] = sorted_id_df[0]

        sample_id_df = pd.read_csv(args.sorted_id_list, delim_whitespace=True, header=None, usecols=[0, 1],
                                   dtype='string')
        sample_id_df["FID"], sample_id_df["IID"] = sample_id_df[0], sample_id_df[1]

        merge_id_df = sorted_id_df.merge(sample_id_df, on="IID", how="inner")
        sample_id_arr = merge_id_df[["FID", "IID"]].to_numpy(dtype="str")
    else:
        #sample_id_arr = np.genfromtxt(sample_file_path, usecols=(0, 1))
        sample_id_arr = np.loadtxt(sample_file_path, usecols=(0, 1), dtype="str")

    if args.randomize:
        print(" Randomizing samples..")
        rng.shuffle(sample_id_arr)

    if args.bisect_ratio:
        write_bisect(sample_id_arr, output_prefix, args.bisect_ratio, n_threads=num_threads)
    elif args.binary_traits:
        if args.overlap_case_samples:
            print(" --overlap-case-samples=true: all control samples will be added to each split")
        #sample_id_df = pd.DataFrame(sample_id_arr, columns=["FID", "IID"], dtype=int)
        sample_id_df = pd.DataFrame(sample_id_arr, columns=["FID", "IID"], dtype='string')
        phen_df = pd.read_csv(args.binary_traits, delim_whitespace=True)
        phen_df[["FID", "IID"]] = phen_df[["FID", "IID"]].astype("string")
        phen_df = sample_id_df.merge(phen_df, on=["FID", "IID"], how="inner")
        cases = phen_df[phen_df.isin([1]).any(axis=1)]
        control_ids = phen_df.index.difference(cases.index)
        print(" Found {}/{} case/control samples of {} total samples across {} traits".format(
            len(cases), len(control_ids), len(phen_df), len(phen_df.columns)-2))
        #control_arr = phen_df.loc[control_ids.to_numpy()][["FID", "IID"]].to_numpy(dtype=int)
        control_arr = phen_df.loc[control_ids.to_numpy()][["FID", "IID"]].to_numpy(dtype="str")
        #case_arr = cases[["FID", "IID"]].to_numpy(dtype=int)
        case_arr = cases[["FID", "IID"]].to_numpy(dtype="str")

        if not args.overlap_case_samples:
            rng.shuffle(case_arr)

        write_splits_bt(control_arr, case_arr, num_splits, output_prefix, overlap_cases=args.overlap_case_samples,
                        n_threads=num_threads)
    else:
        write_splits(sample_id_arr, num_splits, output_prefix, n_threads=num_threads)


if __name__ == "__main__":
    main()