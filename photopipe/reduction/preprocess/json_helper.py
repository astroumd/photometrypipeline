import json


def json_dict_from_file(json_file):
    with open(json_file, 'r') as f:
        json_dict = json.load(f)
    return json_dict


def save_dict_to_json(dictionary, file_name):
    json_string = json.dumps(dictionary, indent=2)
    with open(file_name, 'w') as f:
        f.write(json_string)
