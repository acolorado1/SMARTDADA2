"""
This module focuses contains function wrappers that connects R functions to
python modules.
"""
import glob

from pathlib import Path
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

# import all R modules
_path = Path(__file__)
r_scripts_folder = _path.parent
r_script_path = r_scripts_folder/ "R_scripts"
r_modules: list[str] = glob.glob(f"{r_script_path}/*.R")
r_source = robjects.r["source"]
for r_module in r_modules:
    r_source(r_module)

# # installing dependencies
install_func = robjects.globalenv["RDepInstall"]
install_func()


def r_dist_plot(df) -> None:
    """_summary_

    Parameters
    ----------
    df : _type_
        _description_
    """
    # get R function call
    dist_plot = robjects.globalenv["distribution_boxplox"]

    # call function
    dist_plot(df)


def r_avg_error_plot(df) -> None:
    """_summary_

    Parameters
    ----------
    df : _type_
        _description_
    """

    # convert pandas dataframe in to R dataframe
    r_df = __convert_pandas_to_R_df(df)

    # get R function
    average_error_plot_func = robjects.globalenv["average_error_lineplot"]
    average_error_plot_func(r_df)


def __convert_pandas_to_R_df(df) -> robjects.DataFrame:
    """Conversts pandas dataframe into a R dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        pandas dataframe

    Returns
    -------
    robjects.DataFrame
        R dataframe object
    """

    r_dataframe = pandas2ri.py2rpy(df)
    return r_dataframe