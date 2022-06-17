import sys, time, math
import argparse
import ROOT
import numpy as np
from array import array


def LoadLUT(path_lut):

    print('\nLoadLUT...')

    lut = ROOT.TChain('PhotonLibraryData')
    lut.Add(path_lut)
    lut.SetBranchStatus('*',1)

    print('using lut: ' + path_lut.split('/')[-1] + ' (' + str(lut.GetEntries()) + ' entries)')
    return lut


def LoadDataFile(path_data_file):

    print('\nLoadDataFile...')

    data = ROOT.TChain('t_out')
    data.Add(path_data_file)
    data.SetBranchStatus('*',0)
    data.SetBranchStatus('l_event',1)
    data.SetBranchStatus('l_time',1)
    data.SetBranchStatus('l_tphotons',1)
    data.SetBranchStatus('t_event_ntracks',1)
    data.SetBranchStatus('t_track_start_pos_x',1)
    data.SetBranchStatus('t_track_start_pos_y',1)
    data.SetBranchStatus('t_track_start_pos_z',1)
    data.SetBranchStatus('t_track_end_pos_x',1)
    data.SetBranchStatus('t_track_end_pos_y',1)
    data.SetBranchStatus('t_track_end_pos_z',1)
    data.SetBranchStatus('t_track_length',1)

    print('using data file:  ' + path_data_file.split('/')[-1] + ' (' + str(data.GetEntries()) + ' entries)')
    return data


def GetEntryListData(t_data, event):

    t_data.Draw('>>elist','l_event==%d' % event,'entrylist')

    return(ROOT.gDirectory.Get('elist'))


def GetEntryListLUT(t_lut, voxel):

    t_lut.Draw('>>elist','Voxel==%d' % voxel,'entrylist')

    return(ROOT.gDirectory.Get('elist'))


def GetTpcHeight(path_lut):

    f = ROOT.TFile.Open(path_lut)
    height = f.Max[1]-f.Min[1]
    f.Close()

    return height


def GetTpcDimension(path_lut):

    f = ROOT.TFile.Open(path_lut)
    dim = np.array([f.Max[0]-f.Min[0],f.Max[1]-f.Min[1],f.Max[2]-f.Min[2]])
    f.Close()

    return dim


def GetLutGeometry(path_lut):

    f = ROOT.TFile.Open(path_lut)
    lut_min = np.array([f.Min[0],f.Min[1],f.Min[2]])
    lut_max = np.array([f.Max[0],f.Max[1],f.Max[2]])
    lut_ndiv = np.array([f.NDivisions[0],f.NDivisions[1],f.NDivisions[2]])
    f.Close()

    return np.array([lut_min,lut_max,lut_ndiv])


def GetVoxel(pos,lut_geometry):

    (lut_min,lut_max,lut_ndiv) = lut_geometry
    vox_xyz = np.floor(pos/(lut_max-lut_min)*lut_ndiv).astype(int)+lut_ndiv/2
    voxel = vox_xyz[2]*lut_ndiv[0]*lut_ndiv[1]+vox_xyz[1]*lut_ndiv[0]+vox_xyz[0]

    return voxel


def GetEventNumbers(t_data):

    print('\nGetEventNumbers...')

    event_numbers = []

    for track in t_data:

        if(track.l_event not in event_numbers):
            event_numbers.append(track.l_event)

            sys.stdout.write("\r" + str(len(event_numbers)))
            sys.stdout.flush()

    event_numbers.sort()

    print(' events found!')
    return(event_numbers)


def GetSingleTrackEventNumbers(t_data):

    print('\nGetSingleTrackEventNumbers...')

    event_numbers = []

    for track in t_data:

        if(track.t_event_ntracks!=1):
            continue

        if(track.l_event not in event_numbers):
            event_numbers.append(track.l_event)

            sys.stdout.write("\r" + str(len(event_numbers)))
            sys.stdout.flush()

    event_numbers.sort()

    print(' single-track events found!')
    return(event_numbers)


def GetCosmicTrackEventNumbers(t_data):

    print('\nGetCosmicTrackEventNumbers...')

    event_numbers = []
    track_length = []

    for track in t_data:

        if(track.t_event_ntracks!=1):
            continue
        if(abs(track.t_track_start_pos_y)<145.):
            continue
        if(abs(track.t_track_end_pos_y)<145.):
            continue

        if(track.l_event not in event_numbers):
            event_numbers.append(track.l_event)
            track_length.append(track.t_track_length)

            sys.stdout.write("\r" + str(len(event_numbers)))
            sys.stdout.flush()

    event_numbers.sort()

    print(' cosmic-track candidates found! (average lenght: %.0f mm)' % (sum(track_length)/len(track_length)))
    return(event_numbers)


def GetVerticalCosmicTrackEventNumbers(t_data):

    print('\nGetVerticalCosmicTrackEventNumbers...')

    event_numbers = []
    track_length = []

    for track in t_data:

        if(track.t_event_ntracks!=1):
            continue
        if(abs(track.t_track_start_pos_y)<145.):
            continue
        if(abs(track.t_track_end_pos_y)<145.):
            continue

        track_parameters = GetTrackParameters(track)
        ps = track_parameters[0] # start point
        pe = track_parameters[1] # end point
        pv = track_parameters[2] # pointing vector

        entry_point = ps+300./pv[1]*pv
        exit_point = ps-900./pv[1]*pv

        if(abs(entry_point[0])>150. or abs(exit_point[0])>150.):
            continue
        if(entry_point[2]<0. or exit_point[2]<0. or entry_point[2]>300. or exit_point[2]>300.):
            continue

        if(track.l_event not in event_numbers):
            event_numbers.append(track.l_event)
            track_length.append(track.t_track_length)

            sys.stdout.write("\r" + str(len(event_numbers)))
            sys.stdout.flush()

    event_numbers.sort()

    print(' vertical cosmic-track candidates found! (average lenght: %.0f mm)' % (sum(track_length)/len(track_length)))
    return(event_numbers)


def GetTrackParameters(track):

    start_point = np.array([
        track.t_track_start_pos_x,
        track.t_track_start_pos_y,
        track.t_track_start_pos_z])

    end_point = np.array([
        track.t_track_end_pos_x,
        track.t_track_end_pos_y,
        track.t_track_end_pos_z])

    pointing_vector = end_point-start_point

    pointing_vector = pointing_vector/np.linalg.norm(pointing_vector)

    if(track.t_track_end_pos_y > track.t_track_start_pos_y):
        temp_point = start_point
        start_point = end_point
        end_point = temp_point
        pointing_vector = -1.*pointing_vector

    return([start_point,end_point,pointing_vector])


def LocalToGlobalCoord(vec):
    return(pixelboard_global_coord*np.array([1.,1.,1.])+np.matmul(coord_transform_mask,vec))


def GlobalToLocalCoord(vec):
    return(np.matmul(coord_transform_mask,vec-pixelboard_global_coord*np.array([1.,1.,1.])))


def GetEntryPoint(track_parameters):

    ps_l = track_parameters[0] # start point
    pe_l = track_parameters[1] # end point
    pv_l = track_parameters[2] # pointing vector

    ps_g = LocalToGlobalCoord(ps_l)
    pe_g = LocalToGlobalCoord(pe_l)
    pv_g = np.matmul(coord_transform_mask,pv_l)

    pv_g[pv_g==0] = 1E-9

    tpc_dim = GetTpcDimension(args.lut)
    entry_factor = np.min((ps_g-np.sign(pv_g)*(-1)*tpc_dim/2)/pv_g)
    entry_point = ps_g - entry_factor*pv_g

    return(entry_point)


def GetExitPoint(track_parameters):

    ps_l = track_parameters[0] # start point
    pe_l = track_parameters[1] # end point
    pv_l = track_parameters[2] # pointing vector

    ps_g = LocalToGlobalCoord(ps_l)
    pe_g = LocalToGlobalCoord(pe_l)
    pv_g = np.matmul(coord_transform_mask,pv_l)

    pv_g[pv_g==0] = 1E-9

    tpc_dim = GetTpcDimension(args.lut)
    exit_factor = np.min((np.sign(pv_g)*(1)*tpc_dim/2-pe_g)/pv_g)
    exit_point = pe_g + exit_factor*pv_g

    return(exit_point)

def SimulateEvent(t_lut,t_data,event_number,op_channel_hit_list_data,op_channel_hit_list_lut):

    entry_list_event = GetEntryListData(t_data,event_number)
    print('\nsimulating event no. %d (tree entry %d)...' % (event_number, entry_list_event.GetEntry(0)))

    t_data.GetEntry(entry_list_event.GetEntry(0))

    track_parameters = GetTrackParameters(t_data)

    entry_point = GetEntryPoint(track_parameters)
    exit_point = GetExitPoint(track_parameters)

    print(entry_point)
    print(GlobalToLocalCoord(LocalToGlobalCoord(entry_point)))
    print(LocalToGlobalCoord(entry_point))
    print('\n')
    print(exit_point)
    print(GlobalToLocalCoord(LocalToGlobalCoord(exit_point)))
    print(LocalToGlobalCoord(exit_point))

    track_length = np.linalg.norm(exit_point-entry_point)

    pv = (exit_point-entry_point)/track_length

    lut_geometry = GetLutGeometry(args.lut)

    print('\n  sampling simulated track with step size of %.1f mm...\n' % sampling_step_size)

    voxel = 0
    entry_list_voxel = GetEntryListLUT(t_lut,voxel)

    for step in range(int(track_length/sampling_step_size)+1):
        pos = entry_point+(step+0.5/sampling_step_size)*pv*sampling_step_size

        sys.stdout.write('\r    current position: ' + str(pos) + ('(%.1f %%)' % (np.linalg.norm(entry_point-pos)/track_length*100)))
        sys.stdout.flush()

        voxel_current = GetVoxel(pos,lut_geometry)

        if(voxel!=voxel_current):
            voxel = voxel_current
            entry_list_voxel = GetEntryListLUT(t_lut,voxel)

        for entry in range(entry_list_voxel.GetN()):
            t_lut.GetEntry(entry_list_voxel.GetEntry(entry))
            op_channel_hit_list_lut[t_lut.OpChannel] += t_lut.Visibility*sampling_step_size*photon_yield_mev_mm*energy_deposit_mev_mm

    print('\n\n  simulation summary:')

    for op_channel in sorted(op_channel_map):
        op_channel_hit_list_lut[int(op_channel_map[op_channel])] = round(op_channel_hit_list_lut[int(op_channel_map[op_channel])])
        op_channel_hit_list_data[int(op_channel_map[op_channel])] += round(getattr(t_data, 'l_tphotons/' + op_channel))
        print('    %s: %.0f/%.0f photon hits (lut/true)' % (op_channel,op_channel_hit_list_lut[int(op_channel_map[op_channel])],op_channel_hit_list_data[int(op_channel_map[op_channel])]))

    return 0


def PlotOpResponse(event_numbers,op_channel_hit_list_data,op_channel_hit_list_lut):

    print('\nplotting results...')
    data_file_name = args.data.split('/')[-1].split('.root')[0]
    lut_geometry = GetLutGeometry(args.lut)
    (lut_min,lut_max,lut_ndiv) = lut_geometry

    if(args.n_event>1):
        plot_title = ('op_channel_response')
        plot_name = ('%s_op_channel_response_lut_%d-%d-%d_evt_%d_tot_%d_events' % (data_file_name,lut_ndiv[0],lut_ndiv[1],lut_ndiv[2],event_numbers[0],args.n_event))
    elif(args.event==-999):
        plot_title = ('op_channel_response')
        plot_name = ('%s_op_channel_response_lut_%d-%d-%d_evt_%d' % (data_file_name,lut_ndiv[0],lut_ndiv[1],lut_ndiv[2],event_numbers[0]))
    else:
        plot_title = ('op_channel_response')
        plot_name = ('%s_op_channel_response_lut_%d-%d-%d_evt_%d' % (data_file_name,lut_ndiv[0],lut_ndiv[1],lut_ndiv[2],int(args.event)))

    min_bin = int(op_channel_map[min(op_channel_map, key=op_channel_map.get)])
    max_bin = int(op_channel_map[max(op_channel_map, key=op_channel_map.get)])

    max_hit_value = -999

    hist_lut = ROOT.TH1F(plot_title,'SingleCube-Dec-2020: OpChannel Response;OpChannel [-];Number of Photons [-]',max_bin-min_bin+1,min_bin,max_bin+1)
    hist_data = ROOT.TH1F(plot_title,'SingleCube-Dec-2020: OpChannel Response;OpChannel [-];Number of Photons [-]',max_bin-min_bin+1,min_bin,max_bin+1)

    for op_channel in range(len(op_channel_hit_list_lut)):
        hist_lut.Fill(op_channel,op_channel_hit_list_lut[op_channel])

        if(str(op_channel) in op_channel_map.values() and op_channel_hit_list_lut[op_channel]>max_hit_value):
            max_hit_value = op_channel_hit_list_lut[op_channel]

    for op_channel in range(len(op_channel_hit_list_data)):
        for photon in range(int(op_channel_hit_list_data[op_channel])):
            hist_data.Fill(op_channel)

        if(str(op_channel) in op_channel_map.values() and op_channel_hit_list_data[op_channel]>max_hit_value):
            max_hit_value = op_channel_hit_list_data[op_channel]

    c0 = ROOT.TCanvas('','',1600,1200)
    c0.SetGrid(1)
    c0.SetLeftMargin(0.14)
    c0.SetBottomMargin(0.14)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetLineWidth(4)
    ROOT.gStyle.SetStatX(1.0)
    ROOT.gStyle.SetStatY(0.8)
    ROOT.gStyle.SetStatW(0.2)
    ROOT.gStyle.SetStatH(0.2)
    ROOT.gStyle.SetPalette(103)
    ROOT.gStyle.SetNumberContours(200)
    ROOT.gPad.SetLogz(0)
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()

    hist_lut.SetLineWidth(4)
    hist_lut.SetMaximum(1.1*max_hit_value)
    hist_lut.GetXaxis().SetTitleOffset(1.4)
    hist_lut.GetXaxis().SetTitleSize(0.045)
    hist_lut.GetXaxis().SetLabelSize(0.05)
    hist_lut.GetYaxis().SetTitleOffset(1.5)
    hist_lut.GetYaxis().SetTitleSize(0.045)
    hist_lut.GetYaxis().SetLabelSize(0.05)
    hist_lut.GetZaxis().SetTitleOffset(1.4)
    hist_lut.GetZaxis().SetTitleSize(0.045)
    hist_lut.GetZaxis().SetLabelSize(0.05)
    hist_lut.GetZaxis().RotateTitle(1)

    hist_data.SetLineWidth(4)
    hist_data.SetLineColor(1)
    hist_data.SetMarkerStyle(20)
    hist_data.SetMarkerSize(2)
    hist_data.GetXaxis().SetTitleOffset(1.4)
    hist_data.GetXaxis().SetTitleSize(0.045)
    hist_data.GetXaxis().SetLabelSize(0.05)
    hist_data.GetYaxis().SetTitleOffset(1.5)
    hist_data.GetYaxis().SetTitleSize(0.045)
    hist_data.GetYaxis().SetLabelSize(0.05)
    hist_data.GetZaxis().SetTitleOffset(1.4)
    hist_data.GetZaxis().SetTitleSize(0.045)
    hist_data.GetZaxis().SetLabelSize(0.05)
    hist_data.GetZaxis().RotateTitle(1)

    hist_lut.Draw('hist')
    hist_data.Draw('E1 same')

    legend = ROOT.TLegend(0.70,0.90-2*0.07,0.90, 0.90)
    legend.AddEntry(hist_lut,'lut','l')
    legend.AddEntry(hist_data,'data','lep')
    legend.Draw()

    c0.SaveAs(plot_name+'.png')
    c0.SaveAs(plot_name+'.root')


def PlotSamplingResponse(t_lut,t_data,event_numbers,op_channel_hit_list_data,op_channel_hit_list_lut):

    global sampling_step_size

    sample_op_channel_hit_list_lut = []

    for sampling_step_size in sampling_step_sizes:

        for channel in range(n_op_channel):
            op_channel_hit_list_lut[channel] = 0

        for event_number in event_numbers[0:args.sampling]:
            SimulateEvent(t_lut,t_data,event_number,op_channel_hit_list_data,op_channel_hit_list_lut)

        sample_op_channel_hit_list_lut.append(list(op_channel_hit_list_lut))

    print('\nplotting results...')
    data_file_name = args.data.split('/')[-1].split('.root')[0]
    lut_geometry = GetLutGeometry(args.lut)
    (lut_min,lut_max,lut_ndiv) = lut_geometry

    if(args.sampling==1):
        plot_title = ('op_channel_sampling')
        plot_name = ('%s_op_channel_sampling_lut_%d-%d-%d_evt_%d' % (data_file_name,lut_ndiv[0],lut_ndiv[1],lut_ndiv[2],event_numbers[0]))
    else:
        plot_title = ('op_channel_sampling')
        plot_name = ('%s_op_channel_sampling_lut_%d-%d-%d_evt_%d_tot_%d_events' % (data_file_name,lut_ndiv[0],lut_ndiv[1],lut_ndiv[2],event_numbers[0],args.sampling))

    min_bin = int(op_channel_map[min(op_channel_map, key=op_channel_map.get)])
    max_bin = int(op_channel_map[max(op_channel_map, key=op_channel_map.get)])

    max_hit_value = -999

    hist_list = []

    for sample_list in sample_op_channel_hit_list_lut:

        hist_list.append(ROOT.TH1F(plot_title,'SingleCube-Dec-2020: OpChannel Sampling;OpChannel [-];Number of Photons [-]',max_bin-min_bin+1,min_bin,max_bin+1))

        for op_channel in range(len(sample_list)):
            hist_list[-1].Fill(op_channel,sample_list[op_channel])

            if(str(op_channel) in op_channel_map.values() and sample_list[op_channel]>max_hit_value):
                max_hit_value = sample_list[op_channel]

    c0 = ROOT.TCanvas('','',1600,1200)
    c0.SetGrid(1)
    c0.SetLeftMargin(0.14)
    c0.SetBottomMargin(0.14)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetLineWidth(4)
    ROOT.gStyle.SetStatX(1.0)
    ROOT.gStyle.SetStatY(0.8)
    ROOT.gStyle.SetStatW(0.2)
    ROOT.gStyle.SetStatH(0.2)
    ROOT.gStyle.SetPalette(103)
    ROOT.gStyle.SetNumberContours(200)
    ROOT.gPad.SetLogz(0)
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()

    legend = ROOT.TLegend(0.70,0.90-len(sampling_step_sizes)*0.07,0.90, 0.90)
    hist_key = 0

    for hist in hist_list:
        hist.SetLineWidth(4)
        hist.SetLineColor(ROOT.kRed+4-hist_key)
        hist.SetMaximum(1.1*max_hit_value)
        hist.GetXaxis().SetTitleOffset(1.4)
        hist.GetXaxis().SetTitleSize(0.045)
        hist.GetXaxis().SetLabelSize(0.05)
        hist.GetYaxis().SetTitleOffset(1.5)
        hist.GetYaxis().SetTitleSize(0.045)
        hist.GetYaxis().SetLabelSize(0.05)
        hist.GetZaxis().SetTitleOffset(1.4)
        hist.GetZaxis().SetTitleSize(0.045)
        hist.GetZaxis().SetLabelSize(0.05)
        hist.GetZaxis().RotateTitle(1)

        hist.Draw('hist same')

        legend.AddEntry(hist,'%.0f mm' % sampling_step_sizes[hist_key],'l')

        hist_key += 1

    legend.Draw()

    c0.SaveAs(plot_name+'.png')
    c0.SaveAs(plot_name+'.root')


def WriteTree(t_lut,t_data,event_numbers,op_channel_hit_list_lut,op_channel_t1_list_lut):

    lut_geometry = GetLutGeometry(args.lut)
    tpc_height = GetTpcHeight(args.lut)

    # declare optical channel structure integer
    channel_struct_int = 'struct ChannelListInt {'

    for op_channel in range(n_op_channel):
        channel_struct_int = channel_struct_int + ('int ch%02d; ' % op_channel)

    channel_struct_int = channel_struct_int + '};'

    ROOT.gInterpreter.Declare(channel_struct_int)

    # declare optical channel structure float
    channel_struct_float = 'struct ChannelListFloat {'

    for op_channel in range(n_op_channel):
        channel_struct_float = channel_struct_float + ('float ch%02d; ' % op_channel)

    channel_struct_float = channel_struct_float + '};'

    ROOT.gInterpreter.Declare(channel_struct_float)

    # declare branch variable lists
    l_event = array('i',[0])
    l_time = ROOT.ChannelListFloat()
    l_tphotons = ROOT.ChannelListInt()

    # initialize branch variables
    for op_channel in range(n_op_channel):
        exec('l_time.ch%02d=-999' % op_channel)
        exec('l_tphotons.ch%02d=0' % op_channel)

    # create output tree
    out_file = ROOT.TFile('optSim2x2.root','RECREATE')
    tree = ROOT.TTree('t_out','t_out')

    channel_branch_int = 'ch00/I'

    for op_channel in range(n_op_channel):
        if op_channel == 0:
            continue
        channel_branch_int = channel_branch_int + (':ch%02d' % op_channel)

    channel_branch_float = 'ch00/F'

    for op_channel in range(n_op_channel):
        if op_channel == 0:
            continue
        channel_branch_float = channel_branch_float + (':ch%02d' % op_channel)

    # add branches
    tree.Branch('l_event',l_event,'l_event/I')
    tree.Branch('l_time',l_time,channel_branch_float)
    tree.Branch('l_tphotons',l_tphotons,channel_branch_int)

    # loop events in data file
    for event_number in event_numbers:

        entry_list_event = GetEntryListData(t_data,event_number)
        print('\n-----------------------------------------------------------------------------------------')
        print('simulating event no. %d (tree entry %d), %.1f %%...' % (event_number, entry_list_event.GetEntry(0), float(event_number)/float(event_numbers[-1])*100))

        t_data.GetEntry(entry_list_event.GetEntry(0))

        track_parameters = GetTrackParameters(t_data)

        entry_point = GetEntryPoint(track_parameters)
        exit_point =GetExitPoint(track_parameters)

        track_length = np.linalg.norm(exit_point-entry_point)

        pv = (exit_point-entry_point)/track_length

        lut_geometry = GetLutGeometry(args.lut)

        print('\n  sampling simulated track with step size of %.1f mm...\n' % sampling_step_size)

        voxel = 0
        entry_list_voxel = GetEntryListLUT(t_lut,voxel)

        for op_channel in range(n_op_channel):
            op_channel_hit_list_lut[op_channel] = 0
            op_channel_t1_list_lut[op_channel] = -999

        for step in range(int(track_length/sampling_step_size)+1):
            pos = entry_point+(step+0.5/sampling_step_size)*pv*sampling_step_size

            sys.stdout.write('\r    current position: ' + str(pos) + ('(%.1f %%)' % (np.linalg.norm(entry_point-pos)/track_length*100)))
            sys.stdout.flush()

            voxel_current = GetVoxel(pos,lut_geometry)

            if(voxel!=voxel_current):
                voxel = voxel_current
                entry_list_voxel = GetEntryListLUT(t_lut,voxel)

            for entry in range(entry_list_voxel.GetN()):
                t_lut.GetEntry(entry_list_voxel.GetEntry(entry))

                op_channel = t_lut.OpChannel
                n_photons = t_lut.Visibility*sampling_step_size*photon_yield_mev_mm*energy_deposit_mev_mm
                time = t_lut.T1 + np.linalg.norm(entry_point-pos)/v_c

                op_channel_hit_list_lut[op_channel] += n_photons

                if (op_channel_t1_list_lut[op_channel] == -999 or time < op_channel_t1_list_lut[op_channel]):
                    op_channel_t1_list_lut[op_channel] = time

        for op_channel in range(n_op_channel):
            n_photons = round(op_channel_hit_list_lut[op_channel])
            time = op_channel_t1_list_lut[op_channel]

            exec('l_tphotons.ch%02d = %d' % (op_channel,n_photons))
            exec('l_time.ch%02d = %f' % (op_channel,time))

        l_event[0] = event_number

        tree.Fill()

    tree.Write()


def main():

    t_lut = LoadLUT(args.lut)
    t_data = LoadDataFile(args.data)

    #event_numbers = GetEventNumbers(t_data)
    event_numbers = GetSingleTrackEventNumbers(t_data)
    #event_numbers = GetCosmicTrackEventNumbers(t_data)
    #event_numbers = GetVerticalCosmicTrackEventNumbers(t_data)

    op_channel_hit_list_data = []
    op_channel_hit_list_lut = []
    op_channel_t1_list_lut = []

    for op_channel in range(n_op_channel):
        op_channel_hit_list_data.append(0)
        op_channel_hit_list_lut.append(0)
        op_channel_t1_list_lut.append(-999)

    if(args.list):
        print('\ncosmic track candidate event numbers:')
        print(event_numbers)

    elif(args.tree):

        print('\ncreating optical response root tree...')

        WriteTree(t_lut,t_data,event_numbers,op_channel_hit_list_lut,op_channel_t1_list_lut)

    elif(args.sampling!=-999):

        if(args.sampling>0):
            print('\nsampling %d event(s) at different sampling step-sizes:' % args.sampling)

            if(args.sampling>len(event_numbers)):
                print('\n  only %d suitable events found -> terminating program!' % len(event_numbers))
            else:
                PlotSamplingResponse(t_lut,t_data,event_numbers,op_channel_hit_list_data,op_channel_hit_list_lut)

        else:
            print('\ninvalid sampling number!')

    elif(args.event!=-999):

        if(args.event>=0):
            event_number = args.event

            SimulateEvent(t_lut,t_data,event_number,op_channel_hit_list_data,op_channel_hit_list_lut)

            PlotOpResponse(event_numbers,op_channel_hit_list_data,op_channel_hit_list_lut)

        else:
            print('\ninvalid event number!')

    elif(args.n_event!=-999):

        if(args.n_event>0):
            print('\nsimulating %d event(s):' % args.n_event)

            for event_number in event_numbers[0:args.n_event]:
                SimulateEvent(t_lut,t_data,event_number,op_channel_hit_list_data,op_channel_hit_list_lut)

            PlotOpResponse(event_numbers,op_channel_hit_list_data,op_channel_hit_list_lut)

        else:
            print('\ninvalid numer of events!')

    else:
        event_number = event_numbers[0]
        print('\nno event number was provided --> auto select first cosmic candidate in list (event no. %d)...' % event_number)

        SimulateEvent(t_lut,t_data,event_number,op_channel_hit_list_data,op_channel_hit_list_lut)

        PlotOpResponse(event_numbers,op_channel_hit_list_data,op_channel_hit_list_lut)

    print('\n...no further events to simulate!')


#===============================================>
########## GLOBAL VARIABLE DECLARATION ##########
#===============================================>

v_c =300 # mm/ns

photon_yield_mev_mm = 4000

energy_deposit_mev_mm = 0.2

n_op_channel = 48

op_channel_map = {
        'ch09':'12',
        'ch10':'13',
        'ch11':'14',
        'ch12':'15',
        'ch13':'16',
        'ch14':'17',
        'ch02':'18',
        'ch03':'19',
        'ch04':'20',
        'ch05':'21',
        'ch06':'22',
        'ch07':'23'}

pixelboard_global_coord = np.array([-145.54,155.20,-155.20])    # OpChannel 12-23
coord_transform_mask    = np.array([[0,0,1],[0,1,0],[1,0,0]])   # x- and z-axes flipped

sampling_step_size = 50 # mm

sampling_step_sizes = [300, 100, 30, 10, 3]

#===============================================<


if __name__ == '__main__':
    print('')
    print('===================================================================')
    print('====>  Light Simulation of ArgonCube Cosmic Muon (MIP) Track  <====')
    print('===================================================================')
    print('Patrick Koller, LHEP, University of Bern')
    print('Contact: patrick.koller@lhep.unibe.ch')

    parser = argparse.ArgumentParser(description='Define photon look-up-table, datafile and event to simulate.')
    parser.add_argument('lut', help='path to photon look-up table')
    parser.add_argument('data', help='path to data file')
    parser.add_argument('-e', '--event', type=int, help='number of single event to process', default='-999')
    parser.add_argument('-n', '--n_event', type=int, help='number of events to process', default='-999')
    parser.add_argument('-l', '--list', action='store_true', help='dump cosmic candidate event-numbers')
    parser.add_argument('-t', '--tree', action='store_true', help='create optical response root tree')
    parser.add_argument('-s', '--sampling', type=int, help='compare different sampling step-sizes -> enter number of events to simulate', default='-999')
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(1)

    main()

