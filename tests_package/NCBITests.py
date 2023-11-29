import os
import unittest
from functools import reduce
from hashlib import md5
from inspect import getsourcefile
import operator
from io import StringIO

from Bio import Entrez, SeqIO

from saccharis.NCBIQueries import ncbi_query_dna_from_protein_accessions, download_proteins_from_genomes, \
    ncbi_protein_query

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
test_out_folder = os.path.join(tests_folder, "test_files", "temp")


class NCBITestCase(unittest.TestCase):
    email = "alexscf@msl.ubc.ca"
    accessions = ['EEF89776.1', 'ALJ47680.1', 'EDO56034.1', 'EDO52192.1', 'UBD71357.1', 'UBD71356.1', 'UVP52987.1', 'UVP52988.1', 'QUT88430.1', 'UWZ91296.1', 'ALJ60576.1', 'ALJ60577.1', 'QDM10331.1', 'UWN92831.1', 'UVR10340.1', 'UVO70410.1', 'UVQ30641.1', 'UVR44297.1', 'UVQ38353.1', 'UVR39455.1', 'QUT79263.1', 'QUT79262.1', 'QRQ54585.1', 'QGT72485.1', 'SCV07078.1', 'CAG9870864.1', 'QBJ18101.1', 'QBJ18684.1', 'UVO99996.1', 'UVO99995.1', 'QMI80089.1', 'QMI79558.1', 'UWO00994.1', 'QPH58994.1', 'QNL39874.1', 'UWO08515.1', 'UWO07906.1', 'UVS17843.1', 'UVS18438.1', 'UVS18831.1', 'QUU01445.1', 'QUU00756.1', 'QUT65254.1', 'QUT65987.1', 'QUT61324.1', 'QUT61941.1', 'QUT37093.1', 'QUT35080.1', 'QUT33510.1', 'QQA29804.1', 'QQA30413.1', 'BBK86684.1', 'BBK86024.1', 'UYU51415.1', 'UYU54730.1', 'UYU51622.1', 'QRM99680.1', 'ADB10705.1', 'ADB10706.1', 'QTL79631.1', 'QTL80832.1', 'BAQ28012.1', 'BAQ28011.1', 'QTL78854.1', 'QTL77768.1', 'VEG24685.1', 'VEG24684.1', 'AEY66181.1', 'AEY67147.1', 'CBK96050.1', 'CBK95965.1', 'CBK96026.1', 'CBK96759.1', 'CBK96866.1', 'CBK97722.1', 'UWP26122.1', 'UWP25128.1', 'UWP24391.1', 'UWP26001.1', 'UWP25068.1', 'UWP24506.1', 'UWP24416.1', 'CBL35038.1', 'CBL34359.1', 'CBL35063.1', 'CBL35076.1', 'CBL35246.1', 'CBL35306.1', 'UWP52527.1', 'UWP52744.1', 'UWP52751.1', 'UWP52743.1', 'QNT65311.1', 'QNT65310.1', 'QNT65075.1', 'QNT65320.1', 'UWP57035.1', 'VCV21607.1', 'CBL14187.1', 'QHB23767.1', 'QHB23766.1', 'QEI31268.1', 'QEI31269.1', 'QRT31809.1', 'QRT29985.1', 'AAO75446.1', 'AAC76680.1', 'UBE43496.1', 'UVS32399.1', 'UVS46657.1', 'UVQ72588.1', 'UVR63136.1', 'UVP32767.1', 'UYU58041.1', 'QDM08380.1', 'ALJ44869.1', 'UWN90190.1', 'UBF07300.1', 'UVR07839.1', 'UVO68090.1', 'UVQ33516.1', 'UVQ35485.1', 'UVR36539.1', 'UVQ62767.1', 'UVP10911.1', 'UVP75711.1', 'QUT82294.1', 'WII05610.1', 'QRQ56539.1', 'QGT69668.1', 'SCV09205.1', 'CAG9867637.1', 'UBD66808.1', 'UBD17736.1', 'QUT77841.1', 'UYU43690.1', 'UYU42826.1', 'QBJ18509.1', 'QBJ18351.1', 'QMI79800.1', 'QMI79914.1', 'QPH58649.1', 'QPH58810.1', 'ALJ42575.1', 'UVP60028.1', 'UVV78740.1', 'UVS52656.1', 'UVP54473.1', 'UVR89438.1', 'UVP40611.1', 'UVS11345.1', 'UVV83438.1', 'UVO72197.1', 'UVS23839.1', 'UVP27305.1', 'UVQ40654.1', 'UVQ26412.1', 'UBD10181.1', 'UVQ68111.1', 'QUT70626.1', 'QUT37928.1', 'QZU84228.1', 'QZU78815.1', 'QMW87394.1', 'BCA49648.1', 'QQA06938.1', 'UML59796.1', 'UYV00229.1', 'UYU96782.1', 'UYU92746.1', 'UYU87209.1', 'UYU80595.1', 'UYU78209.1', 'UYU73345.1', 'UYU67439.1', 'UYU60311.1', 'UWO08315.1', 'UWO08178.1', 'UVS18251.1', 'UVS18122.1', 'QUU01228.1', 'QUU01020.1', 'QUT65784.1', 'QUT65612.1', 'QUT61643.1', 'QUT61511.1', 'QUT37350.1', 'QUT37499.1', 'QQA30001.1', 'QQA30149.1', 'BBK86492.1', 'BBK86304.1', 'UYU54558.1', 'QUR43432.1', 'WET86176.1', 'UVQ08905.1', 'UVP22511.1', 'QUT28762.1', 'QUT22975.1', 'QRM97655.1', 'QDH55352.1', 'CBK66673.1', 'AXR42366.1', 'AJE06671.1', 'QGM63544.1', 'BAQ29939.1', 'ADB10715.1', 'ADB10707.1', 'QTL79632.1', 'QTL79640.1', 'BAQ28021.1', 'BAQ28013.1', 'QTL77777.1', 'QTL77769.1', 'VEG24694.1', 'VEG24686.1', 'BAP83969.1', 'AOL10720.1', 'ALE36647.1', 'QNL47153.1', 'QHQ55142.1', 'ACD97519.1', 'QRI56794.1', 'QUI46055.1', 'BCB69086.1', 'QGV03281.1', 'QTL76854.1', 'QWE84212.1', 'QUI44045.1', 'QTL69719.1', 'QWE84287.1', 'WCE37085.1', 'UBS39802.1', 'VEG79810.1', 'BAJ71573.1', 'QOL46912.1', 'QOL29334.1', 'AOP01535.1', 'AOP01536.1', 'UGV93103.1', 'ADQ01987.1', 'UAU38283.1', 'AXF98874.1', 'UNL63124.1', 'ALO75376.1', 'WDY39778.1', 'CBK69798.1', 'QOL41965.1', 'BAJ67056.1', 'QOL38267.1', 'BBV23916.1', 'QSZ16399.1', 'QOL58820.1', 'AEI97107.1', 'QSG89003.1', 'QSY59235.1', 'QOL28577.1', 'UCH54887.1', 'ALO72933.1', 'QLE15601.1', 'UHC28780.1', 'UPW86900.1', 'CBK72749.1', 'UBU24335.1', 'UBU19688.1', 'UBU19885.1', 'UBU21407.1', 'QDK84058.1', 'QZA36518.1', 'QPB32446.1', 'AMG92594.1', 'AVC45122.1', 'BCU46592.1', 'UYF55727.1', 'AKE60660.1', 'AIU00629.1', 'AQS03383.1', 'QUF74946.1', 'UYZ35906.1', 'CUU46017.1', 'QUN36046.1', 'AJG97475.1', 'ABR32956.1', 'ALB47909.1', 'QJU46725.1', 'QCJ08483.1', 'QGH20124.1', 'QGH24159.1', 'QMW92847.1', 'BBK78922.1', 'AXB87220.1', 'WBQ00992.1', 'AEY65035.1', 'AVK49824.1', 'UWP13806.1', 'QRT51213.1', 'AVU18181.1', 'AWX01858.1', 'AWQ41550.1', 'QAZ65103.1', 'QJP75010.1', 'AWQ55767.1', 'ASG40508.1', 'APR43655.1', 'ASB73240.1', 'ASA05478.1', 'AVE72929.1', 'ARA27097.1', 'AVF17599.1', 'AYA14117.1', 'AVO81830.1', 'ASB84954.1', 'AWS80076.1', 'ARZ78878.1', 'AWZ98504.1', 'AKK77239.1', 'AKK90282.1', 'AKK98570.1', 'AKL53941.1', 'QGN40745.1', 'AMY64088.1', 'AKZ71408.1', 'QZS47433.1', 'UAS94702.1', 'ASD57117.1', 'AVP00425.1', 'AWC84073.1', 'AVG36005.1', 'QBN08537.1', 'QUG50264.1', 'QLA68444.1', 'SMF86836.1', 'QLA62367.1', 'AYU96321.1', 'UDG01010.1', 'AIE62021.1', 'AIX52827.1', 'AIX57362.1', 'AEW71524.1', 'QCC94869.1', 'QCC89872.1', 'QCD13198.1', 'QBB07473.1', 'BCT12024.1', 'ASQ75077.1', 'AIV27751.1', 'AMJ68601.1', 'AOE93763.1', 'ASQ15947.1', 'ASO99884.1', 'VEB11776.1', 'AHE72330.1', 'AQT91236.1', 'AVL16644.1', 'QLO45665.1', 'QLU69929.1', 'QFQ07236.1', 'BCZ50701.1', 'BBW29193.1', 'BBW43761.1', 'ADF59672.1', 'AFP68018.1', 'CBK84333.1', 'AFM57873.1', 'AKM85522.1', 'BBS29872.1', 'BBS35211.1', 'BBT42895.1', 'BBT88466.1', 'BBM14795.1', 'AUB54059.1', 'UBM05284.1', 'QCJ55143.1', 'AZP91708.1', 'BAO05609.1', 'SMB21537.1', 'AWY90755.1', 'CAD5438816.1', 'CAD5366009.1', 'CAD5438454.1', 'CAD5450417.1', 'CAD5472561.1', 'CAD6157351.1', 'QBR55266.1', 'CBG36815.1', 'QEG52578.1', 'ANO91605.1', 'ANP09669.1', 'ANP20493.1', 'BCA47000.1', 'QHW61667.1', 'QBR64845.1', 'ARV55469.1', 'AWH68720.1', 'QCN73616.1', 'QCN68163.1', 'API39854.1', 'QCN62686.1', 'QCN57203.1', 'QCN51728.1', 'QCN46254.1', 'ATC15312.1', 'QCN40705.1', 'QCN35154.1', 'QSV49542.1', 'QCN29666.1', 'QCN24180.1', 'QCN18770.1', 'QCN13284.1', 'QCN07796.1', 'QCN02311.1', 'QCM96822.1', 'ATC10416.1', 'QCM91334.1', 'QCM85844.1', 'QCM80361.1', 'ATC09797.1', 'QCM74888.1', 'QCM69414.1', 'QHW72302.1', 'ATB96040.1', 'QCM63924.1', 'AJF58463.1', 'QCM58436.1', 'QCM52949.1', 'QCM47463.1', 'QCM41978.1', 'QCM36491.1', 'QCM31004.1', 'QCM25515.1', 'QCM20026.1', 'QCM14537.1', 'AUY89634.1', 'AUY94790.1', 'AQZ83961.1', 'ARA05422.1', 'ARA00076.1', 'AUY79601.1', 'AUZ08220.1', 'AUZ14805.1', 'AUY85758.1', 'QCL62640.1', 'QCL57152.1', 'QCL51505.1', 'QHW52388.1', 'ATB90905.1', 'QCL45839.1', 'ATB85968.1', 'QCL40269.1', 'QDN14013.1', 'QCL34758.1', 'QCL29325.1', 'QCL23892.1', 'AUK19794.1', 'AUK09365.1', 'AUK04467.1', 'AUJ99008.1', 'AUJ98113.1', 'AUJ89020.1', 'QEZ41569.1', 'QCL18237.1', 'QES62084.1', 'QLK70633.1', 'API00722.1', 'QIL65503.1', 'QHW86086.1', 'QHW81566.1', 'AVU43355.1', 'QHF81792.1', 'AVU63796.1', 'QID73149.1', 'QGL19768.1', 'QGL26694.1', 'QGL31304.1', 'QGL35795.1', 'ATB80999.1', 'QBR14008.1', 'AZH67203.1', 'ALL88237.1', 'ANE59970.1', 'ANE64670.1', 'QHS39580.1', 'AWT03642.1', 'APT64958.1', 'BBP15245.1', 'BBP14941.1', 'BBP22218.1', 'BBP23918.1', 'BBP28349.1', 'BBP32799.1', 'BBP39646.1', 'ARV48954.1', 'BAX13127.1', 'AOR21878.1', 'AVU53694.1', 'QFW20198.1', 'QHF66248.1', 'AVU58540.1', 'AVT61082.1', 'AMF89742.1', 'API06345.1', 'QBY16478.1', 'AMW43442.1', 'QHF71129.1', 'QHI46715.1', 'QDN45320.1', 'ATB76133.1', 'API11879.1', 'QHW43435.1', 'QCX62933.1', 'API17484.1', 'ASZ48121.1', 'QCD71630.1', 'QHW90608.1', 'QCD76121.1', 'AZH61855.1', 'AZH56566.1', 'AZH46071.1', 'AZH41200.1', 'QDN31791.1', 'API23128.1', 'ACI75546.1', 'AZA01428.1', 'AWV13142.1', 'QIR59372.1', 'QHW76735.1', 'ABG71734.1', 'QIZ53382.1', 'QLC53650.1', 'CAV00673.1', 'ARR33007.1', 'QHW57300.1', 'QES63934.1', 'QHW48016.1', 'AJB53678.1', 'ASZ43630.1', 'AVR75293.1', 'AZH35934.1', 'QDM94767.1', 'ATB71075.1', 'API28623.1', 'AJE58225.1', 'QIV72165.1', 'QIV75231.1', 'ART41969.1', 'ACI75547.1', 'ACI75548.1', 'ARV29732.1', 'QHW38850.1', 'API34288.1', 'AIZ89268.1', 'QHI53555.1', 'ART19029.1', 'ART26809.1', 'QDM86179.1', 'AYL87858.1', 'QCV37377.1', 'QCV32310.1', 'QDC28640.1', 'ADN48557.1', 'ASG47544.1', 'AKP86672.1', 'ALY15210.1', 'ATX33617.1', 'QKL27648.1', 'AXF70141.1', 'AXE55336.1', 'AJB36906.1', 'ABJ03135.1', 'AKK35256.1', 'AKK41153.1', 'AGC89142.1', 'QGJ58067.1', 'QEJ88981.1', 'QER56994.1', 'QER61722.1', 'QCR10489.1', 'AWF24511.1', 'AWF20705.1', 'AWF18475.1', 'AWF11426.1', 'ATX41051.1', 'ATX48494.1', 'AXZ41439.1', 'ATX52725.1', 'ATX57293.1', 'ATY20430.1', 'ASC16851.1', 'ARX23745.1', 'AXZ88697.1', 'AVE95949.1', 'AVN37581.1', 'AWS63436.1', 'AXZ01132.1', 'AXZ80793.1', 'ARA66511.1', 'ARZ84002.1', 'ASB78236.1', 'ARW85633.1', 'ARZ87975.1', 'ARX28376.1', 'ARX16522.1', 'AVJ77412.1', 'AVJ70005.1', 'AVN11881.1', 'AWZ82114.1', 'AWZ70207.1', 'AVM04485.1', 'AIL18549.1', 'QJP26086.1', 'QDK08331.1', 'AXV10396.1', 'ACA75741.1', 'AYB95770.1', 'QEF10366.1', 'QEF05723.1', 'QGJ97804.1', 'QGG62087.1', 'AMH24290.1', 'AMH28606.1', 'ACT41157.1', 'AIF63056.1', 'QCT10338.1', 'QDX42237.1', 'QCV93015.1', 'BBU07109.1', 'BBU11534.1', 'CAA0104402.1', 'VZZ90267.1', 'CBV36766.1', 'ATV11021.1', 'AUF93124.1', 'AUQ39542.1', 'AUN92516.1', 'QBO99282.1', 'QBO94571.1', 'QJZ05788.1', 'AJH12285.1', 'CAQ33983.1', 'QJZ14056.1', 'ACT45312.1', 'QHB88099.1', 'ACT27147.1', 'ARD82829.1', 'ARD78961.1', 'ARH99295.1', 'QAU79304.1', 'QAU57880.1', 'QAU62222.1', 'QAU67291.1', 'QAU74579.1', 'QAU80677.1', 'QAU85656.1', 'QAU93456.1', 'QAU95093.1', 'QAV64397.1', 'QAU99387.1', 'QBM95849.1', 'AIN33982.1', 'ACR64142.1', 'ARE45611.1', 'APJ56096.1', 'APJ90111.1', 'APJ96270.1', 'APJ63521.1', 'QJP79079.1', 'QGY26093.1', 'APJ66353.1', 'QOH75717.1', 'APJ70406.1', 'AJF78768.1', 'AKN49522.1', 'AZU43955.1', 'APJ79046.1', 'AXL15749.1', 'APL84614.1', 'APJ84265.1', 'APJ85179.1', 'AWR51519.1', 'QGA34638.1', 'ALX54549.1', 'AST63896.1', 'QHF76376.1', 'AOM55742.1', 'AVS45725.1', 'QBO51320.1', 'QBO40142.1', 'QBO43666.1', 'AKH24324.1', 'ARJ94127.1', 'QEF52426.1', 'QEF56951.1', 'QEF61460.1', 'QEE99387.1', 'QEH87019.1', 'QEH82212.1', 'QEF66079.1', 'QEF70902.1', 'AVS07822.1', 'AAN82925.1', 'NP_756351.1', 'ARQ24140.1', 'AYM21849.1', 'QIB09106.1', 'AKC14192.1', 'QBP03758.1', 'QIF74472.1', 'AUX66187.1', 'AUY75510.1', 'AUY71537.1', 'QKF08126.1', 'AYN88156.1', 'QKF28449.1', 'QKF17782.1', 'QKF12483.1', 'QKF25048.1', 'AXO09591.1', 'ALQ74737.1', 'AUL61476.1', 'AUL67966.1', 'ATX18244.1', 'ATX16990.1', 'ATX10139.1', 'ATW95733.1', 'QFX41779.1', 'AUP44608.1', 'ATZ32732.1', 'ATZ32736.1', 'QEP68156.1', 'QEP73111.1', 'QEP54988.1', 'QEQ81408.1', 'QEQ76450.1', 'QEP59286.1', 'QJH47734.1', 'QEP86838.1', 'QEP63674.1', 'QEP82286.1', 'QEP77869.1', 'APK01549.1', 'APK41813.1', 'APK04355.1', 'APK12951.1', 'APK17591.1', 'APK21343.1', 'APK25131.1', 'APK29998.1', 'QBJ98952.1', 'APK36156.1', 'APK41166.1', 'AWO19455.1', 'AWO27616.1', 'AWO30914.1', 'QJZ10007.1', 'ACX37747.1', 'BAJ45397.1', 'AKR22466.1', 'AKR26819.1', 'AKR31305.1', 'AUT07053.1', 'ASO02759.1', 'AQZ75137.1', 'AYO73590.1', 'AVD30896.1', 'AZZ28759.1', 'BBU77458.1', 'BBF50198.1', 'BBF55856.1', 'BBF62102.1', 'BBU92558.1', 'AXT78887.1', 'AWR71846.1', 'ASL29671.1', 'AXT74347.1', 'QEJ63432.1', 'QCJ58032.1', 'QBA47902.1', 'QBA53352.1', 'AVL08623.1', 'ASA63611.1', 'AZR86992.1', 'QGQ13281.1', 'QAR39607.1', 'QIE79899.1', 'AUO40409.1', 'QDR58706.1', 'QHP97538.1', 'ANO27745.1', 'ASA59652.1', 'AJG10661.1', 'AYP00180.1', 'AZW02542.1', 'AWJ97304.1', 'QDJ85998.1', 'AXZ71305.1', 'QFG19551.1', 'QGQ41509.1', 'ANK50782.1', 'QDM02562.1', 'AQV18366.1', 'AQW19810.1', 'AQV25644.1', 'AQV29797.1', 'AQV35073.1', 'AMX41689.1', 'AQV40224.1', 'AQV46670.1', 'AQV53642.1', 'AQV56976.1', 'AMX13666.1', 'AMX31229.1', 'AMX34238.1', 'AQV63949.1', 'AQV68942.1', 'AQV74121.1', 'AQV81181.1', 'AQV85221.1', 'AQV88128.1', 'AQW00015.1', 'AQW07153.1', 'AQW11868.1', 'AIX65640.1', 'ANK34476.1', 'AUY03391.1', 'AUV31862.1', 'AUV21908.1', 'QRQ31261.1', 'AZV98275.1', 'CAR10475.2', 'AOD09231.1', 'QCW36516.1', 'AYR81362.1', 'ANK07311.1', 'AIZ30145.1', 'QFH28731.1', 'QFH23198.1', 'QFH17922.1', 'QFH12479.1', 'QFH07066.1', 'QFH01554.1', 'QFG96102.1', 'QEO86420.1', 'AXK28329.1', 'AXH83765.1', 'AXH19725.1', 'QFU36899.1', 'CBJ03404.1', 'ATB15823.1', 'ATB10615.1', 'QEP48778.1', 'AUG66831.1', 'QFW01343.1', 'QFW53792.1', 'ARD53454.1', 'AIT36807.1', 'QDB98263.1', 'QQA96474.1', 'QQB96030.1', 'AVG01311.1', 'ATG08348.1', 'ATG12500.1', 'ATM26479.1', 'ATM10505.1', 'ATM82861.1', 'AYZ39832.1', 'AYY94478.1', 'QDD28016.1', 'AUZ91168.1', 'APG34499.1', 'AXI33024.1', 'AOM43112.1', 'APE77663.1', 'AOT30579.1', 'APE89878.1', 'AUX00325.1', 'ASI48366.1', 'ASL56859.1', 'AXV22644.1', 'QBP88269.1', 'QAZ69771.1', 'QGJ25151.1', 'AML06937.1', 'ANO80210.1', 'ASO85572.1', 'QJU28056.1', 'QJU26333.1', 'BBJ89520.1', 'QIB93024.1', 'QEE89674.1', 'QHS06689.1', 'APK48822.1', 'APK81919.1', 'ASA40639.1', 'APK84283.1', 'APK90857.1', 'AWR82017.1', 'ANP34621.1', 'APK55952.1', 'APK61024.1', 'APK64768.1', 'APK71793.1', 'APK74646.1', 'AWJ02935.1', 'AWJ00851.1', 'ARM77701.1', 'QJR23362.1', 'QFV78467.1', 'AVT71138.1', 'QGF36269.1', 'ABV08072.1', 'AUY44809.1', 'AWU88642.1', 'ARM40948.1', 'CCQ31232.2', 'CAR00626.1', 'CAR20293.1', 'QIN55961.1', 'ADE89185.1', 'QGJ10712.1', 'QSS45682.1', 'SCA73545.1', 'SMB21327.1', 'QIF69427.1', 'AVZ47264.1', 'AGY86334.1', 'AMQ53466.1', 'ALX64557.1', 'ALX59768.1', 'QIM34689.1', 'CQR83081.1', 'ALI43054.1', 'QHB67242.1', 'ALI47451.1', 'AMH36981.1', 'AMH32261.1', 'AIZ53469.1', 'AKD93757.1', 'AKD67466.1', 'AKD76185.1', 'AKD84965.1', 'AKD63094.1', 'AKD80596.1', 'AKD89321.1', 'AKD71817.1', 'AOO71863.1', 'AKK14588.1', 'AKK15805.1', 'CDY63007.1', 'CDZ22434.1', 'ANR83479.1', 'QDA48283.1', 'AUY27317.1', 'VWQ03575.1', 'AXE70001.1', 'AIF38942.1', 'ADX48730.1', 'AFH15869.1', 'QCW32318.1', 'ARV34590.1', 'QOH88619.1', 'AZU69866.1', 'QAA00341.1', 'AZQ78352.1', 'QOH83530.1', 'AZU74347.1', 'AZU78834.1', 'AZU85888.1', 'QBG60846.1', 'AYQ10067.1', 'QHL81210.1', 'QHJ62076.1', 'QHI94530.1', 'QHJ54855.1', 'QED50374.1', 'QHJ67957.1', 'CAP78117.1', 'ATV48834.1', 'AGW10710.1', 'QSR37088.1', 'QQX63594.1', 'APK94574.1', 'APL18781.1', 'APL25089.1', 'ATI04704.1', 'APL26289.1', 'ASI14361.1', 'APL33723.1', 'APL39980.1', 'BBG75846.1', 'APK98955.1', 'APL05372.1', 'QFW87969.1', 'APL08210.1', 'AQW75347.1', 'APL15154.1', 'QFX11251.1', 'QMT88589.1', 'CUQ98931.1', 'AXG59139.1', 'AQU94770.1', 'AVZ51795.1', 'ASO76888.1', 'AQP93574.1', 'AUO55282.1', 'AJO85721.1', 'AMB56564.1', 'APE70055.1', 'APE60340.1', 'APE65220.1', 'API49597.1', 'APE55390.1', 'BAX18229.1', 'BAX23105.1', 'AOM72184.1', 'AVQ78312.1', 'AUG95582.1', 'QDJ60062.1', 'QDJ64683.1', 'QAA87883.1', 'AVZ56888.1', 'AML11582.1', 'AEG38648.1', 'AOX55431.1', 'AOX50027.1', 'QIA37280.1', 'QHG54265.1', 'QHG48178.1', 'QHG44654.1', 'QHG39097.1', 'ASQ69290.1', 'AKO54832.1', 'VEC19960.1', 'VDZ29066.1', 'VED38645.1', 'VEC23323.1', 'SQE48476.1', 'VEC14626.1', 'VEA79948.1', 'VEC68372.1', 'VEC36950.1', 'VEC41599.1', 'VEF94773.1', 'VDY73115.1', 'VEC06898.1', 'SNW07259.1', 'VDZ06301.1', 'VEB52586.1', 'VED71685.1', 'AQU01446.1', 'VFQ33423.1', 'VED21865.1', 'CKH01977.1', 'VED00205.1', 'VED04644.1', 'VDY97091.1', 'VDY96477.1', 'VED40925.1', 'VEE27846.1', 'VED09483.1', 'VED09481.1', 'VEC46753.1', 'VEC46752.1', 'VEF19535.1', 'VDY66571.1', 'VEE25910.1', 'VDY91822.1', 'VDY77711.1', 'VDY81956.1', 'VEC76336.1', 'VEC49612.1', 'VED16187.1', 'VEC61991.1', 'VEA49577.1', 'VED23911.1', 'VDY82689.1', 'VEE93218.1', 'VEC64908.1', 'VED31947.1', 'SQE25260.1', 'SQE42084.1', 'QJZ18388.1', 'QJZ26775.1', 'QJY95090.1', 'QJY78955.1', 'QJZ22427.1', 'ANJ40174.1', 'BAU61991.1', 'QJZ01711.1', 'AID80798.1', 'AXY46727.1', 'QAY41417.1', 'QKI62652.1', 'QKI66790.1', 'QKI37801.1', 'AKK44847.1', 'QCJ14451.1', 'QGU41617.1', 'AQX98953.1', 'QFG90840.1', 'QFG85575.1', 'QFG80167.1', 'AWJ34862.1', 'QCH64663.1', 'BAI33202.1', 'ATG64209.1', 'ASF04521.1', 'AVL31663.1', 'QCH54276.1', 'AFS54843.1', 'AFS88731.1', 'AFS72047.1', 'AKE86910.1', 'QKI50918.1', 'AWJ45946.1', 'BBL44047.1', 'BAI38238.1', 'QCH75453.1', 'QCH49220.1', 'QCH44047.1', 'QJE06240.1', 'AWJ27340.1', 'AUF78036.1', 'QCH70109.1', 'SNU19424.1', 'SLM08769.1', 'CAS11476.1', 'ABV19061.1', 'AWN80134.1', 'BBK49640.1', 'BBK54893.1', 'BBK59930.1', 'AHY73139.1', 'AHY67385.1', 'AHG11317.1', 'AHG17062.1', 'QCH91503.1', 'QEQ37831.1', 'ASE48064.1', 'AOV19432.1', 'QKB67661.1', 'QKB78409.1', 'AOV24786.1', 'AOV30137.1', 'QKB94413.1', 'QKB73014.1', 'AOV35504.1', 'AOV40917.1', 'QKB62307.1', 'QKC03687.1', 'QKB89078.1', 'QKB83809.1', 'AOV46262.1', 'QKB56993.1', 'AOV51677.1', 'QKB51698.1', 'QDZ54000.1', 'QKB46238.1', 'QAV75872.1', 'QGG87710.1', 'QDG07118.1', 'QKB41079.1', 'QKB36072.1', 'QKB30950.1', 'QKB26096.1', 'QKB20715.1', 'QKB15373.1', 'QCV05288.1', 'QCV15349.1', 'QKB14620.1', 'QKB09268.1', 'QKA99224.1', 'QKA93850.1', 'QKA88470.1', 'QKA83132.1', 'QKA77848.1', 'QKA72463.1', 'QJZ32748.1', 'QKA61773.1', 'QKA56336.1', 'QKA51135.1', 'ANG82278.1', 'ANG76596.1', 'ANG71100.1', 'QHP64902.1', 'ANW42578.1', 'QCH80468.1', 'QKA45633.1', 'QKA40467.1', 'QKA35125.1', 'QKA29774.1', 'AMG81055.1', 'QKA24389.1', 'QKA18946.1', 'QKA13359.1', 'QKA07853.1', 'QKA02357.1', 'QJZ97108.1', 'QJZ91901.1', 'QJZ86603.1', 'BBC52798.1', 'QJZ81071.1', 'QJZ75668.1', 'QJZ70116.1', 'QJZ64711.1', 'ACI35822.1', 'AAG58801.1', 'NP_290237.1', 'AIG71061.1', 'QKA67085.1', 'BAB37955.1', 'NP_312559.1', 'AIF96222.1', 'AJA28697.1', 'ACT74373.1', 'QJZ59698.1', 'QJZ54361.1', 'AYV43769.1', 'QGF18501.1', 'QJZ48897.1', 'QJZ43504.1', 'ALH92920.1', 'QJZ38109.1', 'QQW38617.1', 'QEJ55637.1', 'QEJ47814.1', 'QEJ43265.1', 'QEJ38954.1', 'QEJ30304.1', 'QEJ34613.1', 'QEJ21606.1', 'QEJ26079.1', 'QEJ19060.1', 'QEJ13227.1', 'QEJ08645.1', 'QEI95422.1', 'QEI89098.1', 'QEI91115.1', 'QEJ04190.1', 'QEI99731.1', 'ANW29648.1', 'QBZ05859.1', 'QDJ69261.1', 'QGQ68085.1', 'ANK04193.1', 'CDN84421.1', 'AWJ57129.1', 'AWJ40982.1', 'QCH97047.1', 'BAI28099.1', 'QEI84266.1', 'AWJ51149.1', 'QCH85842.1', 'ADD58862.1', 'AEZ42764.1', 'AEQ14928.1', 'QEI67313.1', 'ADR29052.1', 'BBM80331.1', 'ATO77633.1', 'QCH59690.1', 'QIL97300.1', 'QIL92125.1', 'AFG42603.1', 'QGY17106.1', 'QGY12373.1', 'QGY08004.1', 'QIM01750.1', 'QGY22095.1', 'APA39915.1', 'ASW61950.1', 'QIF16595.1', 'QIT52017.1', 'QIF11463.1', 'QIB16040.1', 'QIB20413.1', 'ASO90372.1', 'AKK50523.1', 'QHN00010.1', 'QBZ31066.1', 'QCD07182.1', 'QCA20824.1', 'QIQ96321.1', 'QFF93016.1', 'CDH67328.1', 'ASJ45535.1', 'ASJ28600.1', 'ASJ36033.1', 'QGX81667.1', 'QDX34576.1', 'QIE70049.1', 'QIH01684.1', 'QIM39072.1', 'QIE50398.1', 'QLM20020.1', 'QLM24078.1', 'QLM28127.1', 'QLM37321.1', 'QLM41737.1', 'QLM46162.1', 'QLM50889.1', 'QLM93565.1', 'QLN07215.1', 'QLN10490.1', 'QLN12707.1', 'QLN24877.1', 'QLN26073.1', 'QLN31269.1', 'QLO10090.1', 'QLO63683.1', 'QLO66331.1', 'QRB56881.1', 'QRB61702.1', 'QRB38119.1', 'QRA95524.1', 'QRB14486.1', 'QRB03311.1', 'QRA71852.1', 'QRB19008.1', 'QRA81383.1', 'QRB07462.1', 'QRA90920.1', 'QRA86100.1', 'QEG93798.1', 'AWN71180.1', 'QEX65037.1', 'QEX55870.1', 'AVV77521.1', 'AVP31826.1', 'AVV73049.1', 'QGU98552.1', 'QEG88010.1', 'AWZ54197.1', 'QEX68659.1', 'AIZ84696.1', 'AWZ59503.1', 'ALB33677.1', 'AJM75861.1', 'QAU14633.1', 'ALQ57338.1', 'AXA10998.1', 'APL44989.1', 'APL51446.1', 'AXO82769.1', 'QCW48655.1', 'ATZ40305.1', 'AUA44319.1', 'APL57643.1', 'APL49412.1', 'APL61885.1', 'APL64061.1', 'APL72392.1', 'APL75246.1', 'APL79955.1', 'ANJ33086.1', 'APL89320.1', 'QBG44901.1', 'CAR05287.1', 'AZH87457.1', 'AML16604.1', 'AUA41501.1', 'QJI61687.1', 'AYC45221.1', 'AUN45584.1', 'AYL89329.1', 'AUS36206.1', 'QAS83519.1', 'QJP42581.1', 'QJF81880.1', 'QKQ01339.1', 'QJF96120.1', 'QJS68958.1', 'QJG03356.1', 'QJG09895.1', 'QJG10502.1', 'QJG15111.1', 'QJG22114.1', 'QJG24473.1', 'QJG28710.1', 'QJG35625.1', 'QJG38147.1', 'QJG45518.1', 'QJG47609.1', 'QJT88827.1', 'QJF82350.1', 'QJS74377.1', 'QJT83916.1', 'QJG59876.1', 'QJG65075.1', 'QJF91428.1', 'QJS63766.1', 'QJF88443.1', 'QKP92121.1', 'QJB50327.1', 'AWS38934.1', 'BAG79463.1', 'BAI57040.1', 'AKF22834.1', 'ALZ69845.1', 'ALD37904.1', 'ALD33120.1', 'ALD28177.1', 'ALD22946.1', 'QJB63162.1', 'QGJ07113.1', 'BBF21780.1', 'QGJ62413.1', 'APT00552.1', 'AUM20380.1', 'AUL92514.1', 'AUL82850.1', 'ACB16687.1', 'AKF65689.1', 'AKF69829.1', 'AKF73968.1', 'AKF57411.1', 'AKF61551.1', 'AMW48860.1', 'QEM48664.1', 'AHM45907.1', 'AHM50510.1', 'AHM54952.1', 'AXL01648.1', 'AHM27956.1', 'AHM32482.1', 'AHM37044.1', 'ALV71152.1', 'QER73013.1', 'QIP94626.1', 'QHQ99595.1', 'AUM05990.1', 'QHR24650.1', 'QHR29356.1', 'QHR34227.1', 'QHR39304.1', 'QHR44170.1', 'QHR48925.1', 'QHR53461.1', 'QHR58321.1', 'AER91597.1', 'AER86678.1', 'ACB04705.1', 'CDJ73218.1', 'BAL40248.1', 'AMK99849.1', 'APC53756.1', 'BAE77637.1', 'AMU84228.1', 'ARR66141.1', 'QIN68132.1', 'QIN72662.1', 'QJY91059.1', 'QJY82994.1', 'QJY87022.1', 'BCN97405.1', 'ACI75549.1', 'APQ19475.1', 'BCG29823.1', 'BCG34892.1', 'BCG40512.1', 'ACI75550.1', 'QDK73259.1', 'QAY46667.1', 'QAY51359.1', 'QAZ93079.1', 'QBC15335.1', 'ATP24346.1', 'AXP24654.1', 'ALT51540.1', 'ADN73041.1', 'QIH49317.1', 'AEJ59066.1', 'AEE58979.1', 'ANV95706.1', 'AQZ28497.1', 'ABE09637.1', 'AWM67695.1', 'AKA92842.1', 'ADT77272.1', 'AFH13487.1', 'ASO95127.1', 'AVJ15377.1', 'AYL10884.1', 'AVB43577.1', 'AVZ07360.1', 'AYQ03191.1', 'QAS88150.1', 'AVM99178.1', 'QBF89443.1', 'AZR16764.1', 'AWX36508.1', 'APW92805.1', 'AWA14774.1', 'AXN86080.1', 'ATU33026.1', 'BBQ55546.1', 'BBQ60066.1', 'BBQ19159.1', 'BBQ37914.1', 'BBQ42535.1', 'BBQ45947.1', 'BBQ71174.1', 'BBQ92952.1', 'BBQ97823.1', 'BBS65490.1', 'BBS25156.1', 'BBS70412.1', 'BBS70414.1', 'BBS75414.1', 'BBT08918.1', 'BBT23518.1', 'BBT55451.1', 'BBT73617.1', 'BBT97705.1', 'AZM36146.1', 'AZM43409.1', 'QFI48254.1', 'APY02101.1', 'AFJ31326.1', 'APA24298.1', 'ALN47740.1', 'AUO34534.1', 'QEY40770.1', 'QEY36312.1', 'QEY45491.1', 'QEH00294.1', 'ARR58374.1', 'QIG10205.1', 'QIG14476.1', 'AMM38598.1', 'AML21541.1', 'CBK97057.1', 'UWP24807.1', 'CBL33409.1', 'UBD74088.1', 'UPU32821.1', 'APR27609.1', 'QHS02709.1', 'WEE14462.1', 'WIL71758.1', 'WDV25397.1', 'QDJ21710.1', 'AZP89947.1', 'ARW25608.1', 'ARW23607.1', 'WDA28120.1', 'QHM51862.1', 'QHM53885.1', 'QAR86037.1', 'QAT21549.1', 'AOW74873.1', 'UWP57095.1', 'VCV21669.1', 'CBL09656.1', 'UQT32020.1', 'CBL11696.1']
    accessions_sequence_hash = "7972ff99e454dccd8ad216edbbb5f7c8"

    def test_query_dna_from_protein_accessions(self):
        Entrez.email = self.email

        # accessions = ['CCO03766.1', 'CCO03822.1', 'CCO03823.1', 'CCO04221.1', 'CCO04360.1', 'CCO04515.1', 'CCO05195.1', 'CCO05502.1', 'CCO05659.1', 'CCO05987.1', 'CCO06082.1', 'CCO06210.1']
        # accessions = ['CBL18180.1', 'CBL16523.1', 'CBL16847.1', 'CBL16772.1', 'CBL16471.1', 'CBL16630.1', 'CBL16634.1', 'CBL17363.1', 'CBL17440.1', 'CBL17734.1']
        accessions = self.accessions[0:5]
        fasta_dna_records = ncbi_query_dna_from_protein_accessions(accessions)
        fasta_data = reduce(operator.concat, map(lambda record: record.format('fasta'), fasta_dna_records))

        fasta_md5 = md5(fasta_data.encode()).hexdigest()
        self.assertEqual(fasta_md5, self.accessions_sequence_hash)

    def test_protein_query(self):
        fasta_data, queried, retrieved = ncbi_protein_query(self.accessions, api_key=None, ncbi_email=self.email,
                                                            ncbi_tool="saccharis2")
        seqs = list(SeqIO.parse(StringIO(fasta_data), format='fasta'))
        self.assertEqual(len(seqs), len(self.accessions))
        seq_data = ''.join([str(seq.seq) for seq in seqs])
        seq_md5 = md5(seq_data.encode()).hexdigest()
        self.assertEqual(seq_md5, 'cfd595efdb085e0862e83550ab72fd4d')

    def test_query_proteins_from_single_genome(self):
        b_uniformis_genbank = "GCA_018292165.1"
        seqs, sourcedict = download_proteins_from_genomes(b_uniformis_genbank)
        self.assertEqual(len(seqs), 4031)


if __name__ == '__main__':
    unittest.main()
