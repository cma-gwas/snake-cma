from dataclasses import dataclass
import configparser
import os


@dataclass
class ImportConfig:
    """ Class for storing import configuration. """
    name: str
    columns: list[str]
    stat_columns: list[str]
    t_p_cols: list[str]
    is_log10p: bool


def load_config(import_name):
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'config.ini'))
    app_config = config["import.format." + import_name]
    name = import_name
    columns = app_config["columns"].split(',')
    stat_columns = app_config["stat_columns"].split(',')
    t_p_cols = app_config["t_p_columns"].split(',')
    is_log10p = (app_config["is_log10p"].lower() == "true")
    return ImportConfig(name, columns, stat_columns, t_p_cols, is_log10p)

