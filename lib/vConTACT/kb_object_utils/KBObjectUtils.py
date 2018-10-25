import time
import os
import uuid

from KBaseReport.KBaseReportClient import KBaseReport
from KBaseReport.baseclient import ServerError as _RepError
from DataFileUtil.DataFileUtilClient import DataFileUtil as _DFUClient
from DataFileUtil.baseclient import ServerError as _DFUError


def log(message, prefix_newline=False):
    """
    Logging function, provides a hook to suppress or redirect log messages.
    """
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class KBObjectUtils:
    def __init__(self, config):
        self.scratch = os.path.abspath(config['scratch'])
        self.tmp = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(self.tmp)
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.ws_url = config['workspace-url']

    def create_report(self, wsname):
        outdir = os.path.join(self.tmp, 'vcontact_report')
        self._mkdir_p(outdir)

        self._write_search_results(
            os.path.join(outdir, 'index.html'))

        log('Saving Mash search report')

        dfu = _DFUClient(self.callbackURL)
        try:
            dfuout = dfu.file_to_shock({'file_path': outdir, 'make_handle': 0, 'pack': 'zip'})
        except _DFUError as dfue:
            # not really any way to test this block
            log('Logging exception loading results to shock')
            log(str(dfue))
            raise
        log('saved report to shock node ' + dfuout['shock_id'])
        try:
            kbr = KBaseReport(self.callbackURL)
            return kbr.create_extended_report(
                {'direct_html_link_index': 0,
                 'html_links': [{'name': 'index.html',
                                 'label': 'vConTACT results'}
                                ],
                 'report_object_name': 'vcontact_report_' + str(uuid.uuid4()),
                 'workspace_name': wsname
                 })
        except _RepError as re:
            log('Logging exception from creating report object')
            log(str(re))
            # TODO delete shock node
            raise

def _write_search_results(self, outfile):
        # change to mustache or something later. Or just rewrite this whole thing since this is
        # a demo
        with open(outfile, 'w') as html_file:
            html_file.write('<html><body>\n')
            html_file.write('vConTACT Ran successfully')
            html_file.write('</body></html>\n')

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        # https://stackoverflow.com/a/600612/643675
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise
