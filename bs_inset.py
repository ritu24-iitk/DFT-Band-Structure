import numpy as np
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def parse_filband(feig, npl=10):
    # feig : filband in bands.x input file
    # npl : number per line, 10 for bands.x, 6 for phonon

    f=open(feig,'r')
    lines = f.readlines()

    header = lines[0].strip()
    line = header.strip('\n')
    shape = re.split('[,=/]', line)
    nbnd = int(shape[1])
    nks = int(shape[3])
    eig = np.zeros((nks, nbnd), dtype=np.float32)

    dividend = nbnd
    divisor = npl
    div = nbnd // npl + 1 if nbnd % npl == 0 else nbnd // npl + 2
    kinfo=[]
    for index, value in enumerate(lines[1:]):
        value = value.strip(' \n')
        quotient = index // div
        remainder = index % div

        if remainder == 0:
            kinfo.append(value)
        else:
            value = re.split('[ ]+', value)
            a = (remainder - 1) * npl
            b = a + len(value)
            eig[quotient][a:b] = value

    f.close()

    return eig, nbnd, nks, kinfo

def draw_band(bd_file, fig_file, do_find_gap, e_ref=0.0, nvband=0):
    eig, nbnd, nks, kinfo = parse_filband('bands.dat')

    ymin=-1.0  # y range in plot
    ymax=1.0
    lw=1.0 # line width

    p1=plt.subplot(1, 1, 1)
    F=plt.gcf()
    #F.set_size_inches([5,5])

    if nbnd <= nvband:
        print("warning: nvband ", nvband," should be less than the calculated band number ", nbnd)

    plt.xlim([0,nks-1]) # k-points
    plt.ylim([ymin,ymax])
    plt.yticks([-1.0,-0.5,0.0,0.5,1.0], fontsize=20 )
    #plt.xlabel(r'$k (\AA^{-1})$',fontsize=16)
    plt.ylabel(r' E-E$_F$(eV) ',fontsize=20)

    # Set tick marks to be on the inside
    plt.tick_params(axis='both', direction='in')

    if do_find_gap:
        if nbnd > nvband: # for insulators only, nvband can be found by gappw.sh(https://github.com/yyyu200/gappw)
            eig_vbm=max(eig[:,nvband-1])
            eig_cbm=min(eig[:,nvband])
            gap=eig_cbm-eig_vbm
            #plt.title("Band gap= %.4f eV" % (gap))
            e_ref=eig_vbm
        elif nbnd==nvband: # cb not calculated, cannot find gap
            e_ref=max(eig[:,nvband-1])
        else:
            print("set nvband no less than", nbnd)
            assert None
    e_ref=7.4615
    print(e_ref)
    for i in range(nbnd):
        line1=plt.plot( eig[:,i]-e_ref,color='blue',linewidth=lw )
    vlines = [0, 60, 120, 180, 240, 340, 400, 460, 522, 583]

    for vline in vlines:
        plt.axvline(x=vline, ymin=ymin, ymax=ymax, linewidth=lw, color='black')

    xlabeltext = [r'$\Gamma$','L',r'$B_1|B$','Z',r'$\Gamma$',r'$X|Q$','F',r'$P_1$','Z|L','P']

    if len(xlabeltext)<len(vlines):
        for i in range(len(vlines)-len(xlabeltext)):
            xlabeltext.append('X')
    elif len(xlabeltext)>len(vlines):
        xlabeltext=xlabeltext[0:len(vlines)]

    assert len(xlabeltext) == len(vlines)

    plt.xticks( vlines, xlabeltext, fontsize=20 )

    plt.axhline(y=0.0, xmin=0, xmax=nks-1,linewidth=lw,color='black',ls='--' )

    plt.text(4, 8, '$CoNi2Te4$', fontsize=20, color='black', bbox=dict(facecolor='white',alpha=0.99,edgecolor='black') )
    # Increase the left margin slightly more (from 0.15 to 0.2)
    plt.subplots_adjust(left=0.2, bottom=0.15)
    
    # Use tight_layout but tell it to respect the manual adjustment
    plt.tight_layout(rect=[0.05, 0, 1, 1]) 
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

    # --------------------------------------------------
    # Inset around Γ (the second Γ at x = 240)
    # --------------------------------------------------
    ax = plt.gca()

    # position of inset in the main figure
    axins = inset_axes(
        ax,
        width="40%",      # inset width
        height="40%",     # inset height
        bbox_to_anchor=(0.5, 0.04, 1, 1),  # move upward
        bbox_transform=ax.transAxes,
        loc='lower left',
        borderpad=0.5
    )

    # plot all bands again in inset
    for i in range(nbnd):
        axins.plot(
            np.arange(nks),
            eig[:, i] - e_ref,
            color='blue',
            lw=1.0
        )

    # ---- Zoom region ----
    gamma = 240          # Γ position from your vlines list

    axins.set_xlim(gamma-15, gamma+15)
    axins.set_ylim(-0.19, 0.19)

    axins.axvline(
        x=240,
        color='black',
        linewidth=1.0
    )

    axins.axhline(
        y=0.0,
        color='black',
        linestyle='--',
        linewidth=1.0
    )

    # remove ticks for cleaner look
    axins.set_xticks([])
    axins.set_yticks([])

    # inset border thickness
    for spine in axins.spines.values():
        spine.set_linewidth(1.0)
    
    pp, p1, p2 = mark_inset(
    ax,
    axins,
    loc1=1,
    loc2=3,
    fc="none",
    ec="black",
    lw=1.0
    )

    # Dashed connector lines
    p1.set_linestyle('--')
    p2.set_linestyle('--')
    
    plt.savefig(fig_file, dpi=1200)
    

    

if __name__ == '__main__':
    do_find_gap=False
    e_ref=7.4615 # set to fermi-level in scf output for metal, only applicable for do_find_gap=False
    nvband=53 # valence band number, only applicable for do_find_gap=True

    draw_band("bands.dat", "band_BI_inset.png", do_find_gap, nvband=nvband)

