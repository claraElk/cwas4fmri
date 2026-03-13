from cwas4fmri.utils.tools import filter_and_extract_fd
import json
import numpy as np
import tempfile
import os


def test_filter_and_extract_fd():

    sub_list = ["sub-01", "sub-02", "sub-03", "sub-04"]

    # Test the function
    for sub in sub_list:

        # Sub 1 verify that extract FD is correct with 1 run
        if sub == "sub-01":
            # Create a temporary JSON file with mock FD data
            with tempfile.NamedTemporaryFile(
                mode="w+", suffix=".json", delete=False
            ) as tmp_json:
                json_data = {"FDMean": 0.1, "FDMax": 0.2}
                json.dump(json_data, tmp_json)
                tmp_json_path = tmp_json.name

            fdmean = filter_and_extract_fd([tmp_json_path], sub)

            assert fdmean == 0.1  # The mean of the provided FD values

        # Sub 2 verify that subject is excluded based on FD criteria
        elif sub == "sub-02":
            # Create a temporary JSON file
            with tempfile.NamedTemporaryFile(
                mode="w+", suffix=".json", delete=False
            ) as tmp_json:
                json_data = {
                    "FDMean": 0.6,  # Exceeds mean threshold
                    "FDMax": 0.2,
                }
                json.dump(json_data, tmp_json)
                tmp_json_path = tmp_json.name

            fdmean = filter_and_extract_fd([tmp_json_path], sub)

        # Sub 3 verify that subject is excluded based on FD max criteria
        elif sub == "sub-03":
            # Create a temporary JSON file
            with tempfile.NamedTemporaryFile(
                mode="w+", suffix=".json", delete=False
            ) as tmp_json:
                json_data = {
                    "FDMean": 0.1,
                    "FDMax": 3.5,  # Exceeds max threshold
                }
                json.dump(json_data, tmp_json)
                tmp_json_path = tmp_json.name

            fdmean = filter_and_extract_fd([tmp_json_path], sub)

        # Sub 4 verify the mean FD is correctly averaged across multiple runs
        elif sub == "sub-04":
            # Create multiple temporary JSON files with mock FD data
            fd_values = [0.1, 0.2, 0.3]
            json_files = []
            for i, fd in enumerate(fd_values):
                with tempfile.NamedTemporaryFile(
                    mode="w+", suffix=".json", delete=False
                ) as tmp_json:
                    json_data = {"FDMean": fd, "FDMax": 0.2}
                    json.dump(json_data, tmp_json)
                    json_files.append(tmp_json.name)

            fdmean = filter_and_extract_fd(json_files, sub)

            assert fdmean == np.mean(
                fd_values
            )  # The mean of the provided FD values

    # Clean up
    os.remove(tmp_json_path)
