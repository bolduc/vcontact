# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from vConTACT.vConTACTImpl import vConTACT
from vConTACT.vConTACTServer import MethodContext
from vConTACT.authclient import KBaseAuth as _KBaseAuth


class vConTACTTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('vConTACT'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'vConTACT',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = vConTACT(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_vConTACT_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_run_vcontact(self):
        ret = self.getImpl().run_vcontact(self.getContext(), {
            'workspace_name': self.getWsName(),
            # 'genome': '21832/12/1',  # KBaseGenome
            'genome': '21832/10/1',  # KBaseAssembly
            'db': 'ArchaeaViralRefSeq94-Merged',
            'pcs_mode': 'MCL',
            'vcs_mode': 'ClusterONE',
            'blast_evalue': '0.0001',
            'pc_max_overlap': '0.8',
            'pc_penalty': '2',
            'pc_haircut': '0.1',
            'pc_inflation': '2.0',
            'vc_inflation': '2.0',
            'vc_density': '0.3',
            'vc_min_size': '2',
            'vc_max_overlap': '0.9',
            'vc_penalty': '2',
            'vc_haircut': '0.55',
            'merge_method': 'single',
            'similarity': 'match',
            'seed_method': 'nodes',
            'min_significance': '1.0',
            'max_significance': '300',
            'module_inflation': '5.0',
            'mod_significance': '1.0',
            'module_min_shared': '3',
            'link_significance': '1.0',
            'link_proportion': '0.5'
        })
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        pass
