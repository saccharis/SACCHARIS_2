import os
import re
import unittest
from inspect import getsourcefile

from saccharis.ParseUserSequences import parse_multiple_fasta
from saccharis.utils.Formatting import seqs_to_string
from saccharis.utils.PipelineErrors import UserError

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
testfiles_folder = os.path.join(tests_folder, "test_files")
GH102_file = os.path.join(testfiles_folder, "user_test_GH102.fasta")
GH102_UserFormat_file = os.path.join(testfiles_folder, "user_test_GH102_UserFormat.fasta")


def append_to_headers(fasta_data: str, append_string):
    headers = re.findall(">.*\n", fasta_data)
    for header in headers:
        fasta_data = fasta_data.replace(header, f'{header[:-1]}{append_string}\n')
    return fasta_data


class ParseTestCase(unittest.TestCase):

    def test_load_single(self):
        fasta_files = [GH102_file]
        seqs, sources, path = parse_multiple_fasta(fasta_files)
        fasta_string = seqs_to_string(seqs)
        self.assertEqual(fasta_string, FastaData.GH102)

    def test_load_multiple(self):
        fasta_files = [GH102_file, GH102_UserFormat_file]
        seqs, sources, path = parse_multiple_fasta(fasta_files)

        fasta_string = seqs_to_string(seqs)
        fasta_data_string = append_to_headers(FastaData.GH102, f" SACCHARIS merged record from {GH102_file}") + \
                       append_to_headers(FastaData.GH102_UserFormat, f" SACCHARIS merged record from {GH102_UserFormat_file}")
        self.assertEqual(fasta_string, fasta_data_string)

    def test_duplicate_accession(self):
        fasta_files = [GH102_file, GH102_file]
        self.assertRaises(UserError, parse_multiple_fasta, fasta_files)


class FastaData:
    GH102 = ">AAC75420.1 membrane-bound lytic murein transglycosylase A [Escherichia coli str. K-12 substr. MG1655]\nMKGRWVKYLLMGTVVAMLAACSSKPTDRGQQYKDGKFTQPFSLVNQPDAVGAPINAGDFA\nEQINHIRNSSPRLYGNQSNVYNAVQEWLRAGGDTRNMRQFGIDAWQMEGADNYGNVQFTG\nYYTPVIQARHTRQGEFQYPIYRMPPKRGRLPSRAEIYAGALSDKYILAYSNSLMDNFIMD\nVQGSGYIDFGDGSPLNFFSYAGKNGHAYRSIGKVLIDRGEVKKEDMSMQAIRHWGETHSE\nAEVRELLEQNPSFVFFKPQSFAPVKGASAVPLVGRASVASDRSIIPPGTTLLAEVPLLDN\nNGKFNGQYELRLMVALDVGGAIKGQHFDIYQGIGPEAGHRAGWYNHYGRVWVLKTAPGAG\nNVFSGALDVGGAIKGQHFDIY\n>AAW90669.1 murein transglycosylase [Neisseria gonorrhoeae FA 1090]\nMKKHLLRSALYGIAAAILAACQSRSIQTFPQPDTSVINGPDRPAGIPDPAGTTVAGGGAV\nYTVVPHLSMPHWAAQDFAKSLQSFRLGCANLKNRQGWQDVCAQAFQTPVHSFQAKRFFER\nYFTPWQVAGNGSLAGTVTGYYEPVLKGDGRRTERARFPIYGIPDDFISVPLPAGLRGGKN\nLVRIRQTGKNSGTIDNAGGTHTADLSRFPITARTTAIKGRFEGSRFLPYHTRNQINGGAL\nDGKAPILGYAEDPVELFFMHIQGSGRLKTPSGKYIRIGYADKNEHPYVSIGRYMADKGYL\nKLGQTSMQGIKAYMRQNPQRLAEVLGQNPSYIFFRELAGSGNEGPVGALGTPLMGEYAGA\nIDRHYITLGAPLFVATAHPVTRKALNRLIMAQDTGSAIKGAVRVDYFWGYGDEAGELAGK\nQKTTGYVWQLLPNGMKPEYRPWQLLPNGMKPEYRP\n"
    GH102_UserFormat = ">U000000000 AAC75420.1 membrane-bound lytic murein transglycosylase A [Escherichia coli str. K-12 substr. MG1655]\nMKGRWVKYLLMGTVVAMLAACSSKPTDRGQQYKDGKFTQPFSLVNQPDAVGAPINAGDFA\nEQINHIRNSSPRLYGNQSNVYNAVQEWLRAGGDTRNMRQFGIDAWQMEGADNYGNVQFTG\nYYTPVIQARHTRQGEFQYPIYRMPPKRGRLPSRAEIYAGALSDKYILAYSNSLMDNFIMD\nVQGSGYIDFGDGSPLNFFSYAGKNGHAYRSIGKVLIDRGEVKKEDMSMQAIRHWGETHSE\nAEVRELLEQNPSFVFFKPQSFAPVKGASAVPLVGRASVASDRSIIPPGTTLLAEVPLLDN\nNGKFNGQYELRLMVALDVGGAIKGQHFDIYQGIGPEAGHRAGWYNHYGRVWVLKTAPGAG\nNVFSGALDVGGAIKGQHFDIY\n>U000000001 AAW90669.1 murein transglycosylase [Neisseria gonorrhoeae FA 1090]\nMKKHLLRSALYGIAAAILAACQSRSIQTFPQPDTSVINGPDRPAGIPDPAGTTVAGGGAV\nYTVVPHLSMPHWAAQDFAKSLQSFRLGCANLKNRQGWQDVCAQAFQTPVHSFQAKRFFER\nYFTPWQVAGNGSLAGTVTGYYEPVLKGDGRRTERARFPIYGIPDDFISVPLPAGLRGGKN\nLVRIRQTGKNSGTIDNAGGTHTADLSRFPITARTTAIKGRFEGSRFLPYHTRNQINGGAL\nDGKAPILGYAEDPVELFFMHIQGSGRLKTPSGKYIRIGYADKNEHPYVSIGRYMADKGYL\nKLGQTSMQGIKAYMRQNPQRLAEVLGQNPSYIFFRELAGSGNEGPVGALGTPLMGEYAGA\nIDRHYITLGAPLFVATAHPVTRKALNRLIMAQDTGSAIKGAVRVDYFWGYGDEAGELAGK\nQKTTGYVWQLLPNGMKPEYRPWQLLPNGMKPEYRP\n"


if __name__ == '__main__':
    unittest.main()
