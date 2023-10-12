// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnaTaskJpsiVsV0M.h"

// LHC18b
Int_t runList18b[] = { 285447, 285396, 285365, 285364, 285347, 285328, 285327, 285291, 285290, 285289, 285287, 285286, 285224, 285222, 285203, 285202, 285200, 285165, 285127, 285125, 285108, 285106, 285066, 285065, 285064, 285015, 285014, 285013, 285012, 285011, 285010, 285009, 285008 };
// LHC18c
Int_t runList18c[] = { 285958, 285957, 285946, 285917, 285893, 285892, 285869, 285851, 285830, 285812, 285811, 285810, 285806, 285805, 285804, 285781, 285778, 285777, 285756, 285755, 285754, 285753, 285752, 285751, 285722, 285698, 285697, 285664, 285663, 285662, 285659, 285643, 285642, 285641, 285640, 285639, 285603, 285602, 285601, 285599, 285578, 285577, 285576, 285575, 285557, 285515, 285497, 285496 };
// LHC18d
Int_t runList18d[] = { 286350, 286349, 286348, 286345, 286340, 286337, 286336, 286314, 286313, 286312, 286311, 286310, 286309, 286308, 286289, 286288, 286287, 286284, 286282, 286261, 286258, 286257, 286254, 286230, 286229, 286203, 286202, 286201, 286199, 286198, 286159, 286130, 286129, 286127, 286124, 286064, 286028, 286027, 286026, 286025, 286018, 286014, 285980, 285979, 285978 };
// LHC18e
Int_t runList18e[] = { 286937, 286936, 286933, 286932, 286931, 286930, 286911, 286910, 286908, 286907, 286877, 286876, 286874, 286852, 286850, 286848, 286846, 286810, 286809, 286805, 286801, 286799, 286731, 286695, 286661, 286653, 286633, 286594, 286592, 286591, 286569, 286568, 286567, 286566, 286509, 286508, 286502, 286501, 286455, 286454, 286428, 286427, 286426, 286380 };
// LHC18f
Int_t runList18f[] = { 287977, 287975, 287941, 287923, 287784, 287783, 287658, 287657, 287656, 287654, 287578, 287576, 287575, 287573, 287524, 287521, 287520, 287518, 287517, 287516, 287513, 287486, 287484, 287481, 287480, 287451, 287413, 287389, 287388, 287387, 287385, 287381, 287380, 287360, 287358, 287356, 287355, 287353, 287349, 287347, 287346, 287344, 287343, 287325, 287324, 287323, 287283, 287254, 287251, 287250, 287249, 287248, 287209, 287208, 287204, 287203, 287202, 287201, 287155, 287137, 287077, 287072, 287071, 287066, 287064, 287063, 287021, 287000 };
// LHC18l
Int_t runList18l[] = { 289971, 289966, 289943, 289941, 289940, 289935, 289931, 289928, 289888, 289884, 289880, 289857, 289856, 289855, 289852, 289849, 289830, 289816, 289815, 289814, 289811, 289808, 289775, 289757, 289731, 289729, 289724, 289723, 289721, 289666, 289664, 289660, 289659, 289658, 289657, 289654, 289632, 289626, 289625, 289582, 289581, 289579, 289577, 289576, 289574, 289547, 289494, 289493, 289468, 289466, 289465, 289463, 289462, 289444, 289426, 289373, 289370, 289369, 289368, 289367, 289366, 289365, 289363, 289356, 289355, 289354, 289353, 289309, 289308, 289306, 289303, 289300, 289280, 289278, 289277, 289276, 289275, 289254, 289253, 289249, 289247, 289243, 289242, 289241, 289240 };
// LHC18m
Int_t runList18m[] = { 292397, 292298, 292274, 292273, 292270, 292269, 292265, 292242, 292241, 292240, 292192, 292168, 292167, 292166, 292164, 292163, 292162, 292161, 292160, 292140, 292115, 292114, 292109, 292108, 292107, 292106, 292081, 292080, 292077, 292075, 292062, 292061, 292060, 292040, 292012, 291982, 291976, 291953, 291948, 291945, 291944, 291943, 291942, 291803, 291796, 291795, 291769, 291760, 291756, 291755, 291729, 291706, 291698, 291697, 291694, 291692, 291690, 291665, 291661, 291657, 291626, 291625, 291624, 291622, 291618, 291615, 291614, 291590, 291485, 291484, 291482, 291481, 291457, 291456, 291453, 291451, 291447, 291446, 291420, 291419, 291417, 291416, 291402, 291400, 291399, 291397, 291375, 291373, 291363, 291362, 291361, 291360, 291286, 291285, 291284, 291283, 291282, 291265, 291263, 291110, 291100, 291066, 291065, 291041, 291037, 291035, 291006, 291005, 291004, 291003, 291002, 290980, 290979, 290976, 290975, 290948, 290944, 290943, 290935, 290932, 290895, 290894, 290892, 290862, 290860, 290853, 290848, 290790, 290787, 290776, 290774, 290769, 290766, 290764, 290742, 290721, 290699, 290696, 290692, 290687, 290665, 290660, 290658, 290645, 290632, 290627, 290615, 290614, 290613, 290612, 290590, 290553, 290550, 290549, 290544, 290540, 290539, 290538, 290501, 290499, 290469, 290467, 290459, 290458, 290456, 290428, 290427, 290425, 290423, 290421, 290420, 290418, 290411, 290404, 290401, 290375, 290374, 290350, 290327, 290324, 290323, 290300, 290297, 290293, 290254, 290223, 290222 };
// LHC18m - part1
Int_t runList18m1[] = { 292397, 292298, 292274, 292273, 292270, 292269, 292265, 292242, 292241, 292240, 292192, 292168, 292167, 292166, 292164, 292163, 292162, 292161, 292160, 292140, 292115, 292114, 292109, 292108, 292107, 292106, 292081, 292080, 292077, 292075, 292062, 292061, 292060, 292040, 292012, 291982, 291976, 291953, 291948, 291945, 291944, 291943, 291942, 291803, 291796, 291795, 291769, 291760, 291756, 291755, 291729, 291706, 291698, 291697, 291694, 291692, 291690, 291665, 291661, 291657, 291626, 291625, 291624, 291622, 291618, 291615, 291614, 291590, 291485, 291484, 291482, 291481, 291457, 291456, 291453, 291451, 291447, 291446, 291420, 291419, 291417, 291416, 291402, 291400, 291399, 291397, 291375, 291373, 291363, 291362, 291361, 291360 };
// LHC18m - part2
Int_t runList18m2[] = { 291286, 291285, 291284, 291283, 291282, 291265, 291263, 291110, 291100, 291066, 291065, 291041, 291037, 291035, 291006, 291005, 291004, 291003, 291002, 290980, 290979, 290976, 290975, 290948, 290944, 290943, 290935, 290932, 290895, 290894, 290892, 290862, 290860, 290853, 290848, 290790, 290787, 290776, 290774, 290769, 290766, 290764, 290742, 290721, 290699, 290696, 290692, 290687, 290665, 290660, 290658, 290645, 290632, 290627, 290615, 290614, 290613, 290612, 290590, 290553, 290550, 290549, 290544, 290540, 290539, 290538, 290501, 290499, 290469, 290467, 290459, 290458, 290456, 290428, 290427, 290425, 290423, 290421, 290420, 290418, 290411, 290404, 290401, 290375, 290374, 290350, 290327, 290324, 290323, 290300, 290297, 290293, 290254, 290223, 290222 };
// LHC18o
Int_t runList18o[] = { 293898, 293896, 293893, 293891, 293886, 293856, 293831, 293830, 293829, 293809, 293807, 293806, 293805, 293802, 293799, 293776, 293774, 293773, 293741, 293740, 293698, 293696, 293695, 293692, 293691, 293588, 293587, 293497, 293496, 293494, 293475, 293474, 293424, 293413, 293392, 293391, 293388, 293386, 293368 };
// LHC18p
Int_t runList18p[] = { 294925, 294916, 294884, 294883, 294880, 294877, 294875, 294852, 294818, 294817, 294816, 294815, 294813, 294809, 294775, 294774, 294772, 294769, 294749, 294747, 294743, 294742, 294741, 294722, 294721, 294718, 294716, 294715, 294710, 294703, 294653, 294636, 294634, 294633, 294632, 294593, 294591, 294590, 294588, 294587, 294586, 294563, 294558, 294556, 294553, 294531, 294530, 294529, 294527, 294526, 294525, 294524, 294503, 294502, 294310, 294308, 294307, 294305, 294242, 294241, 294212, 294210, 294208, 294205, 294201, 294200, 294199, 294156, 294155, 294154, 294152, 294131, 294128, 294013, 294012, 294011, 294010, 294009 };
// LHC17h
Int_t runList17h[] = { 273103, 273101, 273100, 273099, 273077, 273010, 273009, 272985, 272983, 272976, 272949, 272947, 272939, 272935, 272934, 272933, 272932, 272905, 272903, 272880, 272873, 272871, 272870, 272836, 272835, 272834, 272833, 272829, 272828, 272784, 272783, 272782, 272764, 272763, 272762, 272760, 272755, 272753, 272749, 272747, 272746, 272712, 272692, 272691, 272690, 272620, 272619, 272610, 272608, 272607, 272585, 272577, 272575, 272574, 272521, 272469, 272468, 272466, 272463, 272462, 272461, 272417, 272414, 272413, 272411, 272400, 272399, 272395, 272394, 272389, 272388, 272360, 272359, 272340, 272335, 272194, 272156, 272155, 272154, 272153, 272152, 272151, 272123, 272101, 272100, 272076, 272075, 272042, 272041, 272040, 272039, 272038, 272036, 272034, 272033, 272032, 272031, 272030, 272029, 272025, 272021, 272020, 272018, 271970, 271969, 271962, 271955, 271953, 271946, 271925, 271921, 271916, 271915, 271912, 271911, 271908, 271886, 271881, 271880, 271879, 271878, 271874, 271873, 271871, 271870, 271868 };
// LHC17i
Int_t runList17i[] = { 274442, 274390, 274387, 274385, 274364, 274363, 274360, 274357, 274355, 274329, 274283, 274281, 274280, 274278, 274276, 274271, 274270, 274269, 274268, 274266, 274264, 274263, 274259, 274232, 274212, 274148, 274147, 274125, 274094, 274092, 274064, 274063, 274058, 273986, 273985, 273946, 273942, 273918, 273889, 273887, 273886, 273885, 273825, 273824, 273719, 273711, 273709, 273695, 273690, 273689, 273687, 273654, 273653, 273593, 273592, 273591 };
// LHC17k
Int_t runList17k[] = { 276508, 276507, 276506, 276500, 276462, 276461, 276439, 276438, 276437, 276435, 276434, 276432, 276429, 276351, 276348, 276302, 276297, 276294, 276292, 276291, 276290, 276259, 276230, 276205, 276178, 276177, 276170, 276169, 276166, 276145, 276141, 276140, 276108, 276105, 276104, 276102, 276099, 276098, 275664, 275661, 275657, 275650, 275648, 275624, 275559, 275558, 275515, 275472, 275471, 275467, 275459, 275457, 275453, 275452, 275448, 275406, 275404, 275401, 275369, 275361, 275360, 275357, 275332, 275328, 275283, 275247, 275246, 275245, 275188, 275177, 275175, 275174, 275173, 275151, 275150, 275149, 275076, 275075, 275073, 275070, 275068, 275067, 274979, 274978, 274886, 274884, 274883, 274882, 274822, 274817, 274815, 274811, 274807, 274806, 274803, 274802, 274801, 274743, 274736, 274708 };
// LHC17l
Int_t runList17l[] = { 278216, 278215, 278191, 278189, 278167, 278166, 278165, 278164, 278163, 278130, 278127, 278126, 278123, 278122, 278121, 277996, 277991, 277989, 277988, 277987, 277952, 277930, 277907, 277904, 277903, 277901, 277900, 277899, 277898, 277897, 277876, 277870, 277848, 277847, 277842, 277841, 277836, 277834, 277801, 277800, 277799, 277795, 277794, 277749, 277747, 277746, 277725, 277577, 277576, 277575, 277574, 277537, 277536, 277531, 277530, 277479, 277478, 277476, 277473, 277472, 277470, 277418, 277417, 277389, 277386, 277384, 277383, 277360, 277314, 277312, 277310, 277293, 277262, 277256, 277197, 277196, 277194, 277193, 277189, 277188, 277184, 277183, 277182, 277181, 277180, 277155, 277121, 277117, 277091, 277087, 277082, 277079, 277076, 277073, 277037, 277017, 277016, 277015, 276972, 276971, 276970, 276969, 276920, 276917, 276916, 276762, 276675, 276674, 276672, 276671, 276670, 276669, 276644, 276608, 276557, 276553, 276552, 276551 };
// LHC17m
Int_t runList17m[] = { 280140, 280135, 280134, 280131, 280126, 280118, 280114, 280111, 280108, 280066, 280052, 280051, 280049, 279955, 279954, 279952, 279893, 279890, 279886, 279884, 279880, 279879, 279855, 279854, 279853, 279830, 279827, 279826, 279773, 279749, 279747, 279719, 279718, 279715, 279689, 279688, 279684, 279683, 279682, 279679, 279677, 279676, 279642, 279641, 279600, 279598, 279597, 279583, 279565, 279564, 279563, 279562, 279561, 279560, 279559, 279488, 279487, 279483, 279441, 279439, 279435, 279410, 279391, 279355, 279354, 279349, 279348, 279344, 279342, 279312, 279310, 279309, 279274, 279273, 279270, 279268, 279267, 279265, 279264, 279242, 279238, 279235, 279234, 279208, 279207, 279201, 279199, 279157, 279155, 279130, 279125, 279123, 279122, 279117, 279106, 279075, 279074, 279073, 279068, 279044, 279043, 279041, 279038, 279037, 279036, 279008, 279007, 279005, 278999, 278964, 278963, 278959, 278941, 278939, 278936, 278915, 278914 };
// LHC17o - part1
Int_t runList17o1[] = { 281961, 281956, 281953, 281940, 281939, 281931, 281928, 281918, 281916, 281915, 281894, 281893, 281892, 281755, 281754, 281753, 281751, 281750, 281741, 281713, 281709, 281707, 281706, 281705, 281633, 281592, 281583, 281581, 281580, 281574, 281569, 281568, 281563, 281562, 281557, 281511, 281509, 281477, 281475, 281450, 281449, 281446, 281444, 281441, 281415, 281321, 281301, 281277, 281275, 281244, 281243, 281242, 281241, 281240, 281213, 281212, 281191, 281190, 281181, 281180, 281179, 281081, 281080, 281079, 281062, 281061, 281060, 281036, 281035, 281033, 281032, 280998, 280997, 280996, 280994, 280990, 280947, 280943, 280940, 280936, 280897, 280890, 280881, 280880 };
// LHC17o - part2
Int_t runList17o2[] = { 280856, 280848, 280847, 280845, 280844, 280842, 280793, 280792, 280786, 280768, 280767, 280766, 280765, 280764, 280763, 280761, 280756, 280755, 280754, 280753, 280706, 280705, 280681, 280679, 280676, 280671, 280650, 280648, 280647, 280645, 280639, 280637, 280634, 280613, 280583, 280581, 280576, 280575, 280574, 280551, 280550, 280547, 280546, 280519, 280518, 280448, 280447, 280446, 280445, 280443, 280419, 280415, 280413, 280412, 280406, 280405, 280403, 280375, 280374, 280352, 280351, 280350, 280349, 280348, 280312, 280310, 280290, 280286, 280285, 280284, 280283, 280282 };
// LHC17o - part3
Int_t runList17o3[] = {	281946,281672,281667,281664,281658,281655,281654,281651,281645,281642,281640,281635,281634,280418 };
// LHC17r
Int_t runList17r[] = { 282704, 282703, 282702, 282700, 282677, 282676, 282673, 282671, 282670, 282668, 282667, 282666, 282653, 282651, 282629, 282622, 282620, 282618, 282609, 282608, 282607, 282606, 282580, 282579, 282575, 282573, 282546, 282545, 282544, 282528/*, 282504 */};
// LHC16h
Int_t runList16h[] = { 255469, 255467, 255466, 255465, 255463, 255447, 255442, 255440, 255415, 255402, 255398, 255352, 255351, 255350, 255283, 255280, 255276, 255275, 255256, 255255, 255253, 255252, 255251, 255249, 255248, 255247, 255242, 255240, 255182, 255180, 255177, 255176, 255173, 255171, 255167, 255162, 255159, 255154, 255111, 255091, 255086, 255085, 255082, 255079, 255010, 255009, 255008, 254984, 254983, 254654, 254653, 254652, 254651, 254649, 254648, 254646, 254644, 254640, 254632, 254630, 254629, 254621, 254608, 254606, 254604, 254419 };
// LHC16j
Int_t runList16j[] = { 256420, 256418, 256417, 256415, 256373, 256372, 256371, 256368, 256366, 256365, 256364, 256363, 256362, 256361, 256356, 256311, 256307, 256302, 256298, 256297, 256295, 256292, 256290, 256289, 256287, 256284, 256283, 256282, 256281, 256231, 256228, 256227, 256223, 256222, 256219, 256215, 256213, 256212, 256210, 256204, 256169, 256161, 256158, 256157, 256156, 256149, 256148, 256147, 256146 };
// LHC16k - part1
Int_t runList16k1[] = { 258537, 258499, 258498, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258388, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258280, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258049, 258048, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 257989, 257986, 257979, 257963, 257960, 257958, 257957, 257939, 257937, 257936, 257932, 257912, 257901, 257893, 257892, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257688, 257687, 257685, 257684, 257682, 257644 };
// LHC16k - part2
Int_t runList16k2[] = { 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257565, 257564, 257563, 257562, 257561, 257560, 257541, 257540, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257095, 257092, 257086, 257084, 257083, 257082, 257080, 257077, 257071, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256591, 256567, 256565, 256564, 256561, 256560, 256557, 256556, 256554, 256552, 256512, 256510, 256506, 256504 };
// LHC16o
Int_t runList16o[] = { 264035, 264033, 263985, 263984, 263981, 263979, 263978, 263977, 263923, 263920, 263917, 263916, 263905, 263866, 263863, 263861, 263830, 263829, 263824, 263823, 263813, 263810, 263803, 263793, 263792, 263790, 263787, 263786, 263785, 263784, 263744, 263743, 263741, 263739, 263738, 263737, 263691, 263690, 263689, 263682, 263663, 263662, 263657, 263654, 263653, 263652, 263647, 263529, 263497, 263496, 263490, 263487, 263332, 262858, 262855, 262853, 262849, 262847, 262844, 262842, 262841, 262778, 262777, 262776, 262768, 262760, 262727, 262725, 262723, 262719, 262717, 262713, 262705, 262635, 262632, 262628, 262594, 262593, 262583, 262578, 262574, 262572, 262571, 262570, 262569, 262568, 262567, 262563, 262537, 262533, 262532, 262528, 262492, 262487, 262451, 262430, 262428, 262424, 262423, 262422, 262419, 262418, 262399 };
// LHC16p
Int_t runList16p[] = { 264347, 264346, 264345, 264341, 264336, 264312, 264305, 264281, 264279, 264277, 264273, 264267, 264266, 264265, 264264, 264262, 264261, 264260, 264259, 264238, 264233, 264232, 264198, 264197, 264194, 264188, 264168, 264164, 264139, 264138, 264137, 264129, 264110, 264109, 264086, 264085, 264082, 264078, 264076 };


void runAnalysisONGRID_full(const char *fperiod = "18b", Int_t subset = 0, Bool_t merge = kFALSE, Bool_t mergeJDL = kTRUE, Bool_t local = kFALSE)
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    //Bool_t local = kFALSE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kFALSE;

    //________________Define the periods, the run lists and the path of the data
    TString periodStr = fperiod;
    Int_t nRuns = 0;
    const Int_t *runList;
    if (periodStr.Contains("18b")) { runList = runList18b; nRuns = sizeof(runList18b)/sizeof(runList18b[0]); }
    else if (periodStr.Contains("18c")) { runList = runList18c; nRuns = sizeof(runList18c)/sizeof(runList18c[0]); }
    else if (periodStr.Contains("18d")) { runList = runList18d; nRuns = sizeof(runList18d)/sizeof(runList18d[0]); }
    else if (periodStr.Contains("18e")) { runList = runList18e; nRuns = sizeof(runList18e)/sizeof(runList18e[0]); }
    else if (periodStr.Contains("18f")) { runList = runList18f; nRuns = sizeof(runList18f)/sizeof(runList18f[0]); }
    else if (periodStr.Contains("18l")) { runList = runList18l; nRuns = sizeof(runList18l)/sizeof(runList18l[0]); }
    else if (periodStr.Contains("18m") && subset <= 0) { runList = runList18m; nRuns = sizeof(runList18m)/sizeof(runList18m[0]); }
    else if (periodStr.Contains("18m") && subset == 1) { runList = runList18m1; nRuns = sizeof(runList18m1)/sizeof(runList18m1[0]); }
    else if (periodStr.Contains("18m") && subset == 2) { runList = runList18m2; nRuns = sizeof(runList18m2)/sizeof(runList18m2[0]); }
    else if (periodStr.Contains("18o")) { runList = runList18o; nRuns = sizeof(runList18o)/sizeof(runList18o[0]); }
    else if (periodStr.Contains("18p")) { runList = runList18p; nRuns = sizeof(runList18p)/sizeof(runList18p[0]); }

    else if (periodStr.Contains("17h")) { runList = runList17h; nRuns = sizeof(runList17h)/sizeof(runList17h[0]); }
    else if (periodStr.Contains("17i")) { runList = runList17i; nRuns = sizeof(runList17i)/sizeof(runList17i[0]); }
    else if (periodStr.Contains("17k")) { runList = runList17k; nRuns = sizeof(runList17k)/sizeof(runList17k[0]); }
    else if (periodStr.Contains("17l")) { runList = runList17l; nRuns = sizeof(runList17l)/sizeof(runList17l[0]); }
    else if (periodStr.Contains("17m")) { runList = runList17m; nRuns = sizeof(runList17m)/sizeof(runList17m[0]); }
    else if (periodStr.Contains("17o") && subset == 1) { runList = runList17o1; nRuns = sizeof(runList17o1)/sizeof(runList17o1[0]); }
    else if (periodStr.Contains("17o") && subset == 2) { runList = runList17o2; nRuns = sizeof(runList17o2)/sizeof(runList17o2[0]); }
    else if (periodStr.Contains("17o") && subset == 3) { runList = runList17o3; nRuns = sizeof(runList17o3)/sizeof(runList17o3[0]); }
    else if (periodStr.Contains("17r")) { runList = runList17r; nRuns = sizeof(runList17r)/sizeof(runList17r[0]); }

    else if (periodStr.Contains("16h")) { runList = runList16h; nRuns = sizeof(runList16h)/sizeof(runList16h[0]); }
    else if (periodStr.Contains("16j")) { runList = runList16j; nRuns = sizeof(runList16j)/sizeof(runList16j[0]); }
    else if (periodStr.Contains("16k") && subset == 1) { runList = runList16k1; nRuns = sizeof(runList16k1)/sizeof(runList16k1[0]); }
    else if (periodStr.Contains("16k") && subset == 2) { runList = runList16k2; nRuns = sizeof(runList16k2)/sizeof(runList16k2[0]); }
    else if (periodStr.Contains("16o")) { runList = runList16o; nRuns = sizeof(runList16o)/sizeof(runList16o[0]); }
    else if (periodStr.Contains("16p")) { runList = runList16p; nRuns = sizeof(runList16p)/sizeof(runList16p[0]); }
    else {
      printf("unknown period: %s",fperiod);
      return;
    }

    printf("%d runs will be analyzed in %s\n",nRuns,fperiod);

    TString dataDir;
    if (periodStr.Contains("18")) dataDir = Form("/alice/data/2018/LHC%s",fperiod);
    else if (periodStr.Contains("17")) dataDir = Form("/alice/data/2017/LHC%s",fperiod);
    else if (periodStr.Contains("16")) dataDir = Form("/alice/data/2016/LHC%s",fperiod);
    else {
      printf("unknown year: %s",fperiod);
      return;
    }

    TString dataPattern = "/muon_calo_pass1/AOD/*/AliAOD.root";

    if (periodStr.Contains("17h")) { dataPattern = "/muon_calo_pass2/AOD/*/AliAOD.root"; }
    if (periodStr.Contains("17k")) { dataPattern = "/muon_calo_pass2/AOD/*/AliAOD.root"; }
    if (periodStr.Contains("17l")) { dataPattern = "/pass1/AOD208/*/AliAOD.root"; }
    if (periodStr.Contains("17o") && subset == 1) { dataPattern = "/pass1/AOD208/*/AliAOD.root"; }
    if (periodStr.Contains("17o") && subset == 2) { dataPattern = "/pass1/AOD208/*/AliAOD.root"; }
    if (periodStr.Contains("17o") && subset == 3) { dataPattern = "/muon_calo_pass1/AOD203/*/AliAOD.root"; }
    if (periodStr.Contains("17r")) { dataPattern = "/pass1/AOD208/*/AliAOD.root"; }

    if (periodStr.Contains("16h")) { dataPattern = "/pass1/AOD208/*/AliAOD.root"; }
    if (periodStr.Contains("16j")) { dataPattern = "/pass1/AOD208/*/AliAOD.root"; }
    if (periodStr.Contains("16k")) { dataPattern = "/pass2/AOD208/*/AliAOD.root"; }
    if (periodStr.Contains("16o")) { dataPattern = "/pass1/AOD208/*/AliAOD.root"; }
    if (periodStr.Contains("16p")) { dataPattern = "/pass1/AOD208/*/AliAOD.root"; }

    TString workingDir = Form("LHC%sJpsiV0M_v2",fperiod);
    if (subset > 0) workingDir = Form("LHC%s%dJpsiV0M_v2",fperiod,subset);

    //___________EOF: Define the periods, the run lists and the path of the data


    // since we will compile a class, tell root where to look for headers
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    //loading and executing the task mult selection
    printf("Loading multiplicity selection task...\n");
    TMacro multSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
    AliMultSelectionTask* multSelectionTask = reinterpret_cast<AliMultSelectionTask*>(multSelection.Exec());

    //loading and executing the task physics selection
    printf("Loading physics selection task...\n");
    //TMacro physSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"));
    //AliPhysicsSelectionTask* physSelectionTask = reinterpret_cast<AliPhysicsSelectionTask*>(physSelection.Exec(kFALSE,kTRUE,0,kFALSE));
    //AliPhysicsSelectionTask* physSelectionTask = AddTaskPhysicsSelection(kFALSE,kTRUE,0,kFALSE);
    AliPhysicsSelectionTask *physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(kFALSE,kTRUE,0,kFALSE)"))));

    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnaTaskJpsiVsV0M.cxx++g");
    AliAnaTaskJpsiVsV0M *task = reinterpret_cast<AliAnaTaskJpsiVsV0M*>(gInterpreter->ExecuteMacro("AddJpsiVsV0MTask.C"));
#else
    gROOT->LoadMacro("AliAnaTaskJpsiVsV0M.cxx++g");
    gROOT->LoadMacro("AddJpsiVsV0MTask.C");
    AliAnaTaskJpsiVsV0M *task = AddJpsiVsV0MTask();
#endif


    if(!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);

    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);
    mgr->SetDebugLevel(5);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)
        chain->Add("../myTaskv1/AOD/AliAOD_18b.root");
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {

        if (!TGrid::Connect("alien://")) return;

        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        TString includes_str = "-Wno-deprecated -I$. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include";
        alienHandler->AddIncludePath(includes_str.Data()); // for grid running
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnaTaskJpsiVsV0M.cxx AliAnaTaskJpsiVsV0M.h libPWGmuon.so");
        alienHandler->SetAnalysisSource("AliAnaTaskJpsiVsV0M.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!

        alienHandler->SetMergeViaJDL(mergeJDL);
        if (!merge) {
          alienHandler->SetRunMode("full");
        }
        else {
          alienHandler->SetRunMode("terminate");
        }

        alienHandler->SetAliPhysicsVersion("vAN-20231002_O2-1");//VO_ALICE@AliPhysics::vAN-20231002_O2-1
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data
        alienHandler->SetGridDataDir(dataDir.Data());
        alienHandler->SetDataPattern(dataPattern.Data());

        //TString dataPattern = "/muon_calo_pass1/AOD/*/AliAOD.root";
        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("000");
        // runnumber
        for (Int_t i=0;i<nRuns;i++)  alienHandler->AddRunNumber(runList[i]);
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(75);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(64800);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);

        alienHandler->SetNrunsPerMaster(1);

        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(kFALSE)
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        //alienHandler->SetMergeViaJDL(kTRUE); //do not retrieve final merged output
        //alienHandler->SetMergeViaJDL(kFALSE);  //retrieve final merged output

        // define the output folders
        alienHandler->SetGridWorkingDir(workingDir.Data());
        alienHandler->SetGridOutputDir("Output");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);


        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(21);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            //alienHandler->SetRunMode("full"); //not merging -step1
            //alienHandler->SetRunMode("terminate"); //merging
            mgr->StartAnalysis("grid");
        }
    }
}
