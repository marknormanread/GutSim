"""
Created on 16/04/2014

@author: Mark N. Read
"""

class Nutrients:
    """
    Class provides a convenient way of grouping together nutrient quantities
    """

    def __init__(self, Cw=0., Cd=0., Cs=0., Cm=0., Nf=0., Nm=0., non_ferm=0., Ci=0., other_organic_matter=0.):
        self.Cw = Cw
        self.Cd = Cd
        self.Cs = Cs
        self.Cm = Cm
        self.Nf = Nf
        self.Nm = Nm
        self.Ci = Ci  # Inulin - fibre that Bwm & Bwf can access.
        self.non_fermentable = non_ferm
        # Non-trivial to partition dietary nutrients input to the cecum between their possible destinations. 
        # Some matter will be converted to microbial biomass, and this is broadly retained in the gut.
        # The biomass yield used in Kettle, 2015, Environmental Microbiology 17(5) is around 0.25.
        # The remaining 75% of dietary organic matter is translated into other metabolites.
        # It is unclear how much of this the host will absorb. 
        # This variable captures the portion (in grams) of organic matter that is retained in the gut. 
        # It specifically EXCLUDES the variables accounted for above. 
        self.other_organic_matter = other_organic_matter

    def halve(self):
        """
        Halves the contents of this nutrients object, and returns another Nutrients object with the other half.
        """
        self.Cw /= 2.
        self.Cd /= 2.
        self.Cs /= 2.
        self.Cm /= 2.
        self.Nf /= 2.
        self.Nm /= 2.
        self.Ci /= 2.
        self.non_fermentable /= 2.
        self.other_organic_matter /= 2.
        return Nutrients(Cw=self.Cw, Cd=self.Cd, Cs=self.Cs, Cm=self.Cm, Nf=self.Nf, Nm=self.Nm,
                         Ci=self.Ci, 
                         non_ferm=self.non_fermentable, other_organic_matter=self.other_organic_matter)

    def add(self, other):
        """
        Used for adding various nutrients (as supplied by the `other` Nutrients object) to this object's collection.
        """
        self.Cw += other.Cw
        self.Cd += other.Cd
        self.Cs += other.Cs
        self.Cm += other.Cm
        self.Nf += other.Nf
        self.Nm += other.Nm
        self.Ci += other.Ci
        self.non_fermentable += other.non_fermentable
        self.other_organic_matter += other.other_organic_matter

    def add_quantities(self, Cw=0., Cd=0., Cs=0., Cm=0., Nf=0., Nm=0., non_ferm=0., Ci=0., other_organic_matter=0.):
        self.Cw += Cw
        self.Cd += Cd
        self.Cs += Cs
        self.Cm += Cm
        self.Nf += Nf
        self.Nm += Nm
        self.Ci += Ci
        self.non_fermentable += non_ferm        
        self.other_organic_matter += other_organic_matter

    def subtract(self, other):
        """
        Used for removing various nutrients (as supplied by the `other` Nutrients object) to this object's collection.
        """
        self.Cw -= other.Cw
        self.Cd -= other.Cd
        self.Cs -= other.Cs
        self.Cm -= other.Cm
        self.Nf -= other.Nf
        self.Nm -= other.Nm
        self.Ci -= other.Ci
        self.non_fermentable -= other.non_fermentable
        self.other_organic_matter -= other.other_organic_matter

    def carbon_content(self):
        """ Returns the total quantity of (fermentable) carbon represented in this object. """
        carb = self.Cw + self.Cd + self.Cs + self.Cm + self.Ci
        prot = self.Nf + self.Nm
        C = carbon_content_protein(prot) + carbon_content_carbohydrate(carb)
        return C

    def nitrogen_content(self):
        """ Returns the total quantity of nitrogen represented in this object. """
        prot = self.Nf + self.Nm
        N = nitrogen_content_protein(prot)
        return N

    def total_fermentable_quantity(self):
        """ Returns the total weight of microbe-useable nutrient represented in this object. """
        return self.Cw + self.Cd + self.Cs + self.Cm + self.Nf + self.Nm + self.Ci

    def total_digesta_weight(self):
        return self.total_fermentable_quantity() + self.non_fermentable + self.other_organic_matter

    def prop_Cw(self):
        return self.Cw / self.total_fermentable_quantity()

    def prop_Cd(self):
        return self.Cd / self.total_fermentable_quantity()

    def prop_Cs(self):
        return self.Cs / self.total_fermentable_quantity()

    def prop_Cm(self):
        return self.Cm / self.total_fermentable_quantity()

    def prop_Nf(self):
        return self.Nf / self.total_fermentable_quantity()

    def prop_Nm(self):
        return self.Nm / self.total_fermentable_quantity()

    def prop_inulin(self):
        return self.Ci / self.total_fermentable_quantity()

    def scale_by_factor(self, factor):
        self.Cw *= factor
        self.Cd *= factor
        self.Cs *= factor
        self.Cm *= factor
        self.Nf *= factor
        self.Nm *= factor
        self.Ci *= factor
        self.non_fermentable *= factor
        self.other_organic_matter *= factor

    def proportion_query(self, prop):
        """
        Returns a new Nutrient object containing the proportion of the present object indicated by the argument. 
        The current Nutrient object is unchanged.
        """
        return Nutrients(Cw=prop * self.Cw, 
                         Cd=prop * self.Cd, 
                         Cs=prop * self.Cs, 
                         Cm=prop * self.Cm,
                         Nf=prop * self.Nf, 
                         Nm=prop * self.Nm, 
                         Ci=prop * self.Ci,
                         non_ferm=prop * self.non_fermentable,
                         other_organic_matter=prop * self.other_organic_matter
                         )

    def clone(self):
        """ Returns a copy of this object. """
        return Nutrients(Cw=self.Cw, Cd=self.Cd, Cs=self.Cs, Cm=self.Cm, Nf=self.Nf, Nm=self.Nm,
                         Ci=self.Ci,
                         non_ferm=self.non_fermentable, other_organic_matter=self.other_organic_matter)

    def __str__(self):
        s = "Cw={Cw:.2e} Cd={Cd:.2e} Cm={Cm:.2e} Nf={Nf:.2e} Nm={Nm:.2e} Ci={Ci:.2e} non_ferm={non:.2e} " \
              "oom={oom:.2e} weight={weight:.2e}".format(Cw=self.Cw, Cd=self.Cd, Cm=self.Cm, Nf=self.Nf, Nm=self.Nm,
              Ci=self.Ci, non=self.non_fermentable, oom=self.other_organic_matter, weight=self.total_digesta_weight())
        return(s)


def nitrogen_content_protein(weight):
    """
    Converts the supplied weight of protein into weight of nitrogen.
    """
    # figures are based on these Jones factors:
    # http://www.fao.org/docrep/006/y5022e/y5022e03.htm
    # which suggest that nitrogen comprises around 16% of protein (by mass, or by atoms number? some atoms weigh a lot
    # more than others.).
    # Further support for this figure comes from
    # Rouwenhorst RJ, Jzn JF, Scheffers WA, van Dijken JP. Determination of protein concentration by total organic
    # carbon analysis. Journal of biochemical and biophysical methods. 1991;22(2):119-28.
    return weight / 6.38


def carbon_content_protein(weight):
    """
    Converts the supplied weight of protein into weight of carbon.
    """
    # Evidence for this comes from:
    # Rouwenhorst RJ, Jzn JF, Scheffers WA, van Dijken JP. Determination of protein concentration by total organic
    # carbon analysis. Journal of biochemical and biophysical methods. 1991;22(2):119-28.
    # These authors state that the carbon content of protein is 53% by weight, and that it is stable across different
    # proteins. They do however add that the carbon content of amino acids is 46%, which is different.
    return weight * 0.53


def carbon_content_carbohydrate(weight):
    # evidence for this comes from:
    # Rouwenhorst RJ, Jzn JF, Scheffers WA, van Dijken JP. Determination of protein concentration by total organic
    # carbon analysis. Journal of biochemical and biophysical methods. 1991;22(2):119-28.
    # These authors state that the carbon content of carbohydrate is 44% by weight.
    return weight * 0.44
