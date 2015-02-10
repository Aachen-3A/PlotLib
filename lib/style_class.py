
import matplotlib
import sys

class style_container():
    ##------------------------------------------------------------------
    ## Public functions
    ##------------------------------------------------------------------
    def __init__(self, style = 'Plain', useRoot = False):
        self._style = style

        self._useRoot = useRoot

        #rc('text', usetex=True)
        self._additional_text  = 'Preliminary'
        self._y_label_offset   = -0.11
        self._error_bands_ecol = ['darkmagenta','darkcyan']
        self._error_bands_fcol = ['m','cyan']
        self._error_bands_alph = 0.7
        self._error_bands_labl = ['Sys. uncert. 1','Sys. uncert. 2']
        self._error_bands_center = 'ref'
        self._error_stacking = 'No'
        self._spine_line_width = 0.5
        self._logx = False
        self._logy = True
        self._ymin = -1
        self._ymax = -1
        self._xmin = -1
        self._xmax = -1

    def __del__(self):
        pass

    def AddAxisTitle(self, hist):
        try:
            self._xaxis_title      = hist.xaxis.GetTitle()
            self._yaxis_title      = hist.yaxis.GetTitle()
        except:
            print "Unexpected error:", sys.exc_info()[0], sys.exc_info()[1]
            self._xaxis_title      = 'bla'
            self._yaxis_title      = '#epsilon'

    def InitStyle(self, addplots = ['', '', ''], addheights = [0, 0, 0], cmsPositon = "upper right", legendPosition = "upper right"):
        self.addplots = addplots
        self.addheights = addheights

        self._cmsTextPosition = position(cmsPositon, isText = True)
        self._LegendPosition = position(legendPosition)

        if self._useRoot:
            self._Set_Root_style()
        else:
            matplotlib.rcParams.update({'font.size': 10})
            matplotlib.rcParams.update({'lines.linewidth' : 1})
        if self._style == 'CMS':
            self._Set_CMS_style()
        elif self._style == 'Plain':
            self._Set_Plain_style()
        elif self._style == 'Cool':
            self._Set_Cool_style()
        else:
            print('\n\tNo style chosen, this can\'t be right!\n\tDo something!\n')
            sys.exit(42)

    def Get_useRoot(self):
        return self._useRoot

    def Get_xaxis_title(self):
        return self._xaxis_title

    def Get_yaxis_title(self):
        return self._yaxis_title

    def Get_additional_text(self):
        return self._additional_text

    def Get_y_label_offset(self):
        return self._y_label_offset

    def Get_error_bands_ecol(self):
        return self._error_bands_ecol

    def Get_error_bands_fcol(self):
        return self._error_bands_fcol

    def Get_error_bands_alph(self):
        return self._error_bands_alph

    def Get_error_bands_labl(self):
        return self._error_bands_labl

    def Get_error_bands_center(self):
        return self._error_bands_center

    def Get_error_stacking(self):
        return self._error_stacking

    def Get_spine_line_width(self):
        return self._spine_line_width

    def Get_logx(self):
        return self._logx

    def Get_logy(self):
        return self._logy

    def Get_ymin(self):
        return self._ymin

    def Get_ymax(self):
        return self._ymax

    def Get_xmin(self):
        return self._xmin

    def Get_xmax(self):
        return self._xmax

    def Get_add_cms_text(self):
        return self._add_cms_text

    def Get_add_lumi_text(self):
        return self._add_lumi_text

    def Get_label_text_color(self):
        return self._label_text_color

    def Get_annotation_text_color(self):
        return self._annotation_text_color

    def Get_bg_color(self):
        return self._bg_color

    def Get_ref_line_color(self):
        return self._ref_line_color

    def Get_spine_color(self):
        return self._spine_color

    def Get_tick_color(self):
        return self._tick_color

    def Get_marker_style(self):
        return self._marker_style

    def Get_marker_size(self):
        return self._marker_size

    def Get_marker_color(self):
        return self._marker_color

    def Get_marker_error_cap_width(self):
        return self._marker_error_cap_width

    def Get_cms_text_alignment(self):
        return self._cms_text_alignment

    def Get_show_minor_tick_labels(self):
        return self._show_minor_tick_labels

    def Get_legend_font_size(self):
        return self._legend_font_size

    def Get_cmsTextPosition(self):
        return self._cmsTextPosition

    def Get_LegendPosition(self):
        return self._LegendPosition



    def Set_error_bands_labl(self, label):
        self._error_bands_labl = label

    def Set_error_bands_center(self, center):
        self._error_bands_center = center

    def Set_error_stacking(self, stacking):
        self._error_stacking = stacking

    ## Function to set properties of the plotting axis
    #
    # This function sets axis properties like the y-range or
    # if any axis should be logarithmic.
    # @param[in] logx Boolean if the x-axis should be logarithmic (Default = False)
    # @param[in] logy Boolean if the y-axis should be logarithmic (Default = True)
    # @param[in] ymin Minimum plotting range for the y-axis (Default = -1 automatic values)
    # @param[in] ymax Maximum plotting range for the y-axis (Default = -1 automatic values)
    # @param[in] xmin Minimum plotting range for the x-axis (Default = -1 range from hist)
    # @param[in] xmax Maximum plotting range for the x-axis (Default = -1 range from hist)
    def Set_axis(self, logx = False, logy = True, ymin = -1, ymax = -1, xmin = -1, xmax = -1):
        self._logx = logx
        self._logy = logy
        self._ymin = ymin
        self._ymax = ymax
        self._xmin = xmin
        self._xmax = xmax
    ##------------------------------------------------------------------
    ## Private functions
    ##------------------------------------------------------------------
    def _Set_CMS_style(self):
        self._add_cms_text           = True
        self._add_lumi_text          = True
        self._label_text_color       = 'black'
        self._annotation_text_color  = 'black'
        self._bg_color               = 'w'
        self._ref_line_color         = 'blue'
        self._spine_color            = 'black'
        self._tick_color             = 'black'
        self._marker_style           = 'o'
        self._marker_size            = 3
        self._marker_color           = 'black'
        self._marker_error_cap_width = 0
        self._cms_text_alignment     = 'row'
        self._show_minor_tick_labels = False
        self._legend_font_size       = 9
        if self.addplots[0] != '':
            self._cmsTextPosition.addYspace(  -0.9 * self.addheights[0] / 100.)
        if self.addplots[1] != '':
            self._cmsTextPosition.addYspace(  0.9 * self.addheights[1] / 100.)
        if self.addplots[2] != '':
            self._cmsTextPosition.addYspace(  0.9 * self.addheights[2] / 100.)

    def _Set_Plain_style(self):
        self._add_cms_text           = False
        self._add_lumi_text          = False
        self._label_text_color       = 'black'
        self._annotation_text_color  = 'black'
        self._bg_color               = 'w'
        self._ref_line_color         = 'blue'
        self._spine_color            = 'black'
        self._tick_color             = 'black'
        self._marker_style           = 'o'
        self._marker_size            = 4
        self._marker_color           = 'black'
        self._marker_error_cap_width = 1
        self._cms_text_alignment     = 'row'
        self._show_minor_tick_labels = True
        self._legend_font_size       = 10
        if self.addplots[0] != '':
            self._cmsTextPosition.addYspace(  -0.8 * self.addheights[0] / 100.)
        if self.addplots[1] != '':
            self._cmsTextPosition.addYspace(  0.8 * self.addheights[1] / 100.)
        if self.addplots[2] != '':
            self._cmsTextPosition.addYspace(  0.8 * self.addheights[2] / 100.)

    def _Set_Cool_style(self):
        self._add_cms_text           = True
        self._add_lumi_text          = True
        self._label_text_color       = 'white'
        self._annotation_text_color  = 'white'
        self._bg_color               = '#07000d'
        self._ref_line_color         = 'y'
        self._spine_color            = '#5998ff'
        self._tick_color             = 'w'
        self._marker_style           = 'o'
        self._marker_size            = 3
        self._marker_color           = 'lightgray'
        self._marker_error_cap_width = 0
        self._cms_text_alignment     = 'column'
        self._show_minor_tick_labels = False
        self._legend_font_size       = 9

    def _Set_Root_style(self):
        self.cmsTextFont          = 61   # Fonts
        self.lumiTextFont         = 42
        self.extraTextFont        = 52
        self.additionalTextFont   = 42
        self.cmsTextSize          = 0.9  #Text sizes
        self.lumiTextSize         = 0.6
        self.extraTextSize        = 0.76*self.cmsTextSize
        self.additionalTextSize   = 1.0*self.extraTextSize
        self.legendTextSize       = self.extraTextSize*0.8
        self.lumiTextOffset       = 0.2
        self.extraTextOffset      = 2.5  # only used in outOfFrame version
        self.axisTextSize         = 0.9
        self.axisOffset           = 1.3
        self._ratio_pad           ={}
        self.rootMemory           =[]

class position():
    def __init__(self,positiontext="upper right", refference="", isText=False):

        self._positiontext=positiontext
        if not isinstance(positiontext,str):
            self._definedCoorinates=True
            self._valign="left"
            self._align="left"
        else:
            self._definedCoorinates=False
            self._valign=self._positiontext.split(" ")[0]
            self._align=self._positiontext.split(" ")[1]
        self.addY=0
        self.addX=0
        self._isText=isText
        self._correctcms={"left":0.,
                    "middle":0.,
                    "right":-0.15,
                    "upper":-0.04,
                    "center":0.,
                    "lower":0.,
        }
        if self._definedCoorinates:
            self._x=self._positiontext[0]
            self._y=self._positiontext[1]

    def __eq__(self,other):
        return (self._positiontext==other._positiontext)

    def addYspace(self,y):
        if self._valign=="upper" and y<0.:
            self.addY+=y
        elif self._valign!="upper" and self.getY()+y>0.1:
            self.addY+=y

    def addXspace(self,x):
        self.addX+=x



    def setPosition(self,positiontext):
        self._positiontext=positiontext
        self.valign=self._positiontext.split(" ")[0]
        self.align=self._positiontext.split(" ")[1]

    def getText(self):
        return self._positiontext

    def getX(self):
        if self._definedCoorinates:
            return self._x
        alignDict={
                    "left":0.12,
                    "middle":0.5,
                    "right":0.95,
        }
        if self._isText:
            return self.addX+alignDict[self._align]+self._correctcms[self._align]
        return self.addX+alignDict[self._align]

    def getY(self):
        if self._definedCoorinates:
            return self._y
        alignDict={
                    "upper":0.95,
                    "center":0.5,
                    "lower":0.12,
        }
        if self._isText:
            return self.addY+alignDict[self._valign]+self._correctcms[self._valign]
        return self.addY+alignDict[self._valign]
